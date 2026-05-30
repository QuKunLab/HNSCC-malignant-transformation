#!/usr/bin/env python3
import subprocess
from pathlib import Path
import shlex
import sys
import json
import xml.etree.ElementTree as ET
import unicodedata
import tifffile as tiff



CORES_OME_DIR = Path("/cores_ome")  
SEG_ROOT = Path("/segmentation_results")  
OUT_BASE = Path("/h5ad_out") 
QupathMeasurement2h5ad_SCRIPT = Path("/src/QupathMeasurement2h5ad.py")


RESOLUTION = 0.325
SITE = "Cell"
STAT = "Mean"
NO_CIRCLE_CLIP = False 

RETURN_CELL_MASK = True
RETURN_BOUNDARY_MASK = True

python_exe = sys.executable


def get_hw(ome_path: Path):
    with tiff.TiffFile(str(ome_path)) as tf:
        shp = tf.series[0].shape
    if len(shp) == 3:
        if shp[0] <= 200:  # (C,Y,X)
            return int(shp[1]), int(shp[2])
        else:              # (Y,X,C)
            return int(shp[0]), int(shp[1])
    elif len(shp) == 2:
        return int(shp[0]), int(shp[1])
    raise ValueError(f"Unexpected shape for {ome_path.name}: {shp}")


def norm_name(s: str) -> str:
    s = unicodedata.normalize("NFKC", s)
    s = s.replace(" ", "").replace("\t", "").strip().lower()
    return s


def read_ome_channel_names(ome_path: Path):
    with tiff.TiffFile(str(ome_path)) as tf:
        ome_xml = tf.ome_metadata
    if not ome_xml:
        return []
    root = ET.fromstring(ome_xml)
    ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
    names = []
    for ch in root.findall(".//ome:Channel", ns):
        nm = ch.get("Name") or ""
        if nm:
            names.append(nm)
    return names


def read_geojson_channels(geojson_path: Path, site=SITE, stat=STAT):
    with open(geojson_path, "r") as f:
        data = json.load(f)
    feats = [x for x in data.get("features", []) if x.get("properties", {}).get("objectType") == "cell"]
    if not feats:
        return []
    meas = feats[0].get("properties", {}).get("measurements", {})
    suffix = f": {site}: {stat}"
    chans = []
    for k in meas.keys():
        if k.endswith(suffix):
            chans.append(k[: -len(suffix)])
    return chans


def infer_channels_ome_order_but_geojson_names(ome_path: Path, geojson_path: Path, site=SITE, stat=STAT):
    ome_ch = read_ome_channel_names(ome_path)
    gj_ch = read_geojson_channels(geojson_path, site=site, stat=stat)
    if not gj_ch:
        return []

    gj_map = {}
    for c in gj_ch:
        k = norm_name(c)
        if k not in gj_map:
            gj_map[k] = c  

    selected = []
    used = set()
    for c in ome_ch:
        k = norm_name(c)
        if k in gj_map and gj_map[k] not in used:
            selected.append(gj_map[k])
            used.add(gj_map[k])

    return selected if selected else gj_ch


def core_id_from_ome(ome_path: Path) -> str:
    """1-A.ome.tif -> 1-A"""
    name = ome_path.name
    for suf in [".ome.tif", ".ome.tiff", ".tif", ".tiff"]:
        if name.endswith(suf):
            return name[:-len(suf)]
    return ome_path.stem


def geojson_path_from_core_id(core_id: str) -> Path:
    p = SEG_ROOT / f"{core_id}.ome.tif - Image0" / "cell_seg" / "detections.geojson"
    if not p.exists():
        raise FileNotFoundError(f"Missing geojson: {p}")
    return p


def main():
    OUT_BASE.mkdir(parents=True, exist_ok=True)

    ome_files = sorted(CORES_OME_DIR.glob("*.ome.tif*"))
    if not ome_files:
        raise SystemExit(f"No *.ome.tif(f) found in {CORES_OME_DIR}")

    n_ok = n_fail = n_skip = 0

    for ome in ome_files:
        core_id = core_id_from_ome(ome)
        core_out = OUT_BASE / core_id
        core_out.mkdir(parents=True, exist_ok=True)

        out_h5ad = core_out / f"{core_id}.h5ad"
        if out_h5ad.exists():
            print(f"✅ Skip (exists): {out_h5ad}")
            n_skip += 1
            continue

        try:
            geojson = geojson_path_from_core_id(core_id)
            H, W = get_hw(ome)

            channels = infer_channels_ome_order_but_geojson_names(ome, geojson, site=SITE, stat=STAT)
            if not channels:
                raise RuntimeError(f"[{core_id}] Cannot infer channels (no measurements keys found?).")

            cmd = [
                python_exe, str(QupathMeasurement2h5ad_SCRIPT),
                "-j", str(geojson),
                "-f", str(out_h5ad),
                "-s", str(H), str(W),
                "--offset_px", "0", "0",        
                "-r", str(RESOLUTION),
                "--site", SITE, "--stat", STAT,
                "--return_cell_mask",
                "--return_boundary_mask",
                "-c", *channels,
            ]

            if NO_CIRCLE_CLIP:
                cmd.append("--no_circle_clip")

            print(f"\n[{core_id}] HxW={H}x{W}, channels={len(channels)} (first5={channels[:5]})")
            print("Running:", " ".join(shlex.quote(x) for x in cmd))
            subprocess.run(cmd, check=True)
            n_ok += 1

        except Exception as e:
            print(f"❌ Failed core={core_id}: {e}")
            n_fail += 1

    print(f"\nALL DONE. ok={n_ok}, fail={n_fail}, skip={n_skip}")


if __name__ == "__main__":
    main()
