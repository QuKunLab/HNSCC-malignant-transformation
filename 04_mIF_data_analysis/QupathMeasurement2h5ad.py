#!/usr/bin/env python3
import argparse
import os
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import cv2
import anndata as ad
from scipy.sparse import csr_matrix
from tifffile import imwrite

warnings.filterwarnings("ignore")



def parse_args():
    parser = argparse.ArgumentParser(description="GeoJSON (QuPath) to cell mask + h5ad (no detections.tsv)")

    parser.add_argument("-j", "--geojson", type=str, required=True,
                        help="Full path of QuPath exported GEOJSON file (per-core / WSI).")

    parser.add_argument("-f", "--h5ad_path", type=str, default=f"{os.getcwd()}/adata.h5ad",
                        help="Output h5ad path. Default: ./adata.h5ad")

    parser.add_argument("-s", "--image_shape", type=int, nargs=2, required=True,
                        metavar=("H", "W"),
                        help="Image shape (H W) in pixels for generating masks (core ROI / WSI).")

    parser.add_argument("-c", "--channel", type=str, nargs="+", required=True,
                        help="Protein channels to keep (names must match GEOJSON measurements keys).")

    parser.add_argument("-r", "--resolution", type=float, default=0.325,
                        help="Resolution: um per pixel. Default 0.325")

    parser.add_argument("--offset_px", type=int, nargs=2, default=[0, 0],
                        metavar=("X1", "Y1"),
                        help="Offset in FULL-image pixel coords: (x1 y1). Used to shift polygons into core ROI.")

    parser.add_argument("--site", type=str, default="Cell",
                        choices=["Cell", "Nucleus", "Cytoplasm", "Membrane"],
                        help="Which compartment to use in measurements. Default: Cell")
    parser.add_argument("--stat", type=str, default="Mean",
                        choices=["Mean", "Median", "Min", "Max", "Std.Dev."],
                        help="Which statistic to use in measurements. Default: Mean")

    parser.add_argument("--core_diameter_mm", type=float, default=2.4,
                        help="Core diameter in mm. Default 2.4")
    parser.add_argument("--no_circle_clip", action="store_true",
                        help="Disable circular clipping for masks and cells.")

    parser.add_argument("--return_cell_mask", action="store_true",
                        help="Save cell_mask.tif next to h5ad.")
    parser.add_argument("--return_boundary_mask", action="store_true",
                        help="Save boundary_mask.tif next to h5ad.")

    parser.add_argument("--thickness", type=int, default=1,
                        help="Boundary thickness (pixels). Default 1")

    return parser.parse_args()



def polygon_centroid_xy(points_xy: np.ndarray) -> tuple[float, float]:
    """
    points_xy: (N,2), polygon vertices in pixel coords.
    """
    x = points_xy[:, 0]
    y = points_xy[:, 1]
    # ensure closed
    if x[0] != x[-1] or y[0] != y[-1]:
        x = np.r_[x, x[0]]
        y = np.r_[y, y[0]]

    a = x[:-1] * y[1:] - x[1:] * y[:-1]
    A = a.sum() / 2.0
    if abs(A) < 1e-6:
        return float(points_xy[:, 0].mean()), float(points_xy[:, 1].mean())

    cx = ((x[:-1] + x[1:]) * a).sum() / (6.0 * A)
    cy = ((y[:-1] + y[1:]) * a).sum() / (6.0 * A)
    return float(cx), float(cy)



def generate_masks_from_geojson(
        geojson_file: str,
        *,
        image_size: tuple[int, int],  # (H, W) core ROI
        offset_px: tuple[int, int] = (0, 0),  # (x1,y1) FULL-image px -> local
        thickness: int = 1,
        enable_circle_clip: bool = True,
        circle_center: tuple[int, int] | None = None,  # (cx, cy) local px
        circle_radius_px: int | None = None,
        return_cell_id: bool = True,
):
    H, W = image_size
    cell_mask = np.zeros((H, W), dtype=np.int32)
    boundary_mask = np.zeros((H, W), dtype=np.uint8)

    with open(geojson_file, "r") as f:
        data = json.load(f)

    cells = [x for x in data["features"] if x.get("properties", {}).get("objectType") == "cell"]

    ox, oy = offset_px

    if circle_center is None:
        circle_center = (W // 2, H // 2)
    if circle_radius_px is None:
        circle_radius_px = min(H, W) // 2

    cx0, cy0 = circle_center
    r2 = circle_radius_px * circle_radius_px

    kept = []
    label = 1

    for feat in cells:
        obj_id = feat.get("id")
        geom = feat.get("geometry", {})

        if geom.get("type") == "Polygon":
            coords = geom["coordinates"][0]
        else:
            coords = geom["coordinates"][0][0]  # MultiPolygon: first polygon

        pts = np.asarray(coords, dtype=np.float32)
        # shift to local ROI
        pts[:, 0] -= ox
        pts[:, 1] -= oy

        if enable_circle_clip:
            mx, my = pts[:, 0].mean(), pts[:, 1].mean()
            if (mx - cx0) ** 2 + (my - cy0) ** 2 > r2:
                continue

        # clip to ROI
        pts[:, 0] = np.clip(pts[:, 0], 0, W - 1)
        pts[:, 1] = np.clip(pts[:, 1], 0, H - 1)
        pts_i32 = pts.astype(np.int32).reshape((-1, 1, 2))

        cv2.fillPoly(cell_mask, [pts_i32], label)
        cv2.polylines(boundary_mask, [pts_i32], isClosed=True, color=255, thickness=thickness)

        if return_cell_id:
            kept.append((label, obj_id))

        label += 1

    if enable_circle_clip:
        yy, xx = np.ogrid[:H, :W]
        circle = (xx - cx0) ** 2 + (yy - cy0) ** 2 <= r2
        cell_mask[~circle] = 0
        boundary_mask[~circle] = 0

    if return_cell_id:
        df_cell_id = pd.DataFrame(kept, columns=["cell_id", "Object ID"])
        return cell_mask, boundary_mask, df_cell_id
    return cell_mask, boundary_mask


def geojson_measurements_to_h5ad(
        geojson_file: str,
        *,
        channel_list: list[str],
        site: str,
        stat: str,
        resolution_um_per_px: float,
        offset_px: tuple[int, int],
        image_size: tuple[int, int],
        enable_circle_clip: bool,
        circle_radius_px: int,
):
    H, W = image_size
    ox, oy = offset_px

    with open(geojson_file, "r") as f:
        data = json.load(f)

    feats = [x for x in data["features"] if x.get("properties", {}).get("objectType") == "cell"]

    expr_keys = [f"{ch}: {site}: {stat}" for ch in channel_list]

    X_rows = []
    obs_rows = []
    spatial_rows = []
    obj_ids = []

    # circle filter params (local)
    cx0, cy0 = W // 2, H // 2
    r2 = circle_radius_px * circle_radius_px

    morph_keys = [
        "Detection probability",
        "Nucleus: Area µm^2", "Nucleus: Length µm", "Nucleus: Circularity", "Nucleus: Solidity",
        "Nucleus: Max diameter µm", "Nucleus: Min diameter µm",
        "Cell: Area µm^2", "Cell: Length µm", "Cell: Circularity", "Cell: Solidity",
        "Cell: Max diameter µm", "Cell: Min diameter µm",
        "Nucleus/Cell area ratio",
    ]

    for feat in feats:
        obj_id = feat.get("id")
        meas = feat.get("properties", {}).get("measurements", {})


        if any(k not in meas for k in expr_keys):
            continue

        # polygon -> local coords
        geom = feat.get("geometry", {})
        if geom.get("type") == "Polygon":
            coords = geom["coordinates"][0]
        else:
            coords = geom["coordinates"][0][0]

        pts = np.asarray(coords, dtype=np.float32)
        pts[:, 0] -= ox
        pts[:, 1] -= oy

        cx, cy = polygon_centroid_xy(pts)


        if not (0 <= cx < W and 0 <= cy < H):
            continue

        if enable_circle_clip:
            if (cx - cx0) ** 2 + (cy - cy0) ** 2 > r2:
                continue

        expr = [meas[k] for k in expr_keys]
        X_rows.append(expr)

        # obs
        obs = {
            "Object ID": obj_id,
            "pixel_X": int(round(cx)),
            "pixel_Y": int(round(cy)),
            "Centroid X µm": float(cx * resolution_um_per_px),
            "Centroid Y µm": float(cy * resolution_um_per_px),
        }
        for k in morph_keys:
            if k in meas:
                obs[k] = meas[k]
        obs["signal_sum"] = float(np.sum(expr))
        obs_rows.append(obs)

        spatial_rows.append([int(round(cx)), int(round(cy))])
        obj_ids.append(obj_id)

    if len(obj_ids) == 0:
        raise ValueError(
            "No cells parsed from geojson. "
            "Check channel names + site/stat keys exist (e.g. 'HIF-1α: Cell: Mean'), and offset/image_shape."
        )

    X = csr_matrix(np.asarray(X_rows, dtype=np.float32))
    adata = ad.AnnData(X=X)
    adata.var_names = pd.Index(channel_list, name="protein")
    adata.obs = pd.DataFrame(obs_rows, index=pd.Index(obj_ids, name="cell_id"))
    adata.obsm["spatial"] = np.asarray(spatial_rows, dtype=np.float32)


    adata.uns["expr_source"] = {"site": site, "stat": stat}
    adata.uns["resolution_um_per_px"] = float(resolution_um_per_px)
    adata.uns["offset_px"] = {"x": int(ox), "y": int(oy)}
    adata.uns["image_shape_hw"] = {"H": int(H), "W": int(W)}

    return adata


def main():
    args = parse_args()

    H, W = args.image_shape
    RESOLUTION = float(args.resolution)  # um/px
    offset_px = (int(args.offset_px[0]), int(args.offset_px[1]))

    out_h5ad = Path(args.h5ad_path)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)

    # core circle radius (px)
    core_radius_um = (float(args.core_diameter_mm) * 1000.0) / 2.0
    circle_radius_px = int(round(core_radius_um / RESOLUTION))
    enable_circle_clip = (not args.no_circle_clip)

    # 1) masks from geojson polygons
    cell_mask, boundary_mask, df_cell_id = generate_masks_from_geojson(
        args.geojson,
        image_size=(H, W),
        offset_px=offset_px,
        thickness=int(args.thickness),
        enable_circle_clip=enable_circle_clip,
        circle_center=(W // 2, H // 2),
        circle_radius_px=circle_radius_px,
        return_cell_id=True,
    )

    if args.return_cell_mask:
        imwrite(str(out_h5ad.parent / "cell_mask.tif"), cell_mask.astype(np.int32), compression="zlib")
    if args.return_boundary_mask:
        imwrite(str(out_h5ad.parent / "boundary_mask.tif"), boundary_mask.astype(np.uint8), compression="zlib")

    # 2) h5ad from geojson measurements (no tsv)
    adata = geojson_measurements_to_h5ad(
        args.geojson,
        channel_list=args.channel,
        site=args.site,
        stat=args.stat,
        resolution_um_per_px=RESOLUTION,
        offset_px=offset_px,
        image_size=(H, W),
        enable_circle_clip=enable_circle_clip,
        circle_radius_px=circle_radius_px,
    )

    oid2cid = dict(zip(df_cell_id["Object ID"], df_cell_id["cell_id"]))

    adata.obs["cell_id_label"] = adata.obs["Object ID"].map(oid2cid)

    valid_mask = adata.obs["cell_id_label"].notna()
    n_dropped = (~valid_mask).sum()

    if n_dropped > 0:
        print(f"[Warning] Dropping {n_dropped} cells that exist in measurements but not in mask (boundary mismatch).")
        adata = adata[valid_mask].copy()

    adata.obs["cell_id_label"] = adata.obs["cell_id_label"].astype(int)
    adata.obs_names = pd.Index([f"cell_{x}" for x in adata.obs["cell_id_label"]], name="cell")

    adata.uns["core_diameter_mm"] = float(args.core_diameter_mm)
    adata.uns["core_radius_um"] = float(core_radius_um)
    adata.uns["circle_radius_px"] = int(circle_radius_px)
    adata.uns["circle_clip_enabled"] = bool(enable_circle_clip)

    adata.write_h5ad(str(out_h5ad))


if __name__ == "__main__":
    main()