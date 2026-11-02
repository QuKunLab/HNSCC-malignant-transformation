library(dplyr)
library(survival)
library(plyr)
library(dplyr)
library(survival)
library(plyr)
library(tools) 
cox_pipe <- function(data,
                     var_n,
                     filter_by_p_value = FALSE,  
                     p_cutoff = 0.05,           
                     covariate_order = NULL      
){
  
  colnames(data)[2] <- "status"
  colnames(data)[3] <- "time"
  
  surv <- Surv(time = data$time, event = data$status)
  data$surv <- surv

  default_vars <- colnames(data)[4:var_n]
  
  if (is.null(covariate_order)) {
    varNames <- default_vars
  } else {
    if (!all(covariate_order %in% default_vars)) {
      stop("covariate_order 中有不在数据协变量列里的名字，请检查。")
    }
    varNames <- covariate_order
  }

  UniCox <- function(x) {
    FML <- as.formula(paste0("surv ~ ", x))
    Cox <- coxph(FML, data = data)
    Sum <- summary(Cox)
    
    CI <- paste0(round(Sum$conf.int[, 3:4], 2), collapse = "-")
    Pvalue <- round(Sum$coefficients[, 5], 3)
    HR <- round(Sum$coefficients[, 2], 2)
    
    Unicox <- data.frame(
      "Characteristics" = x,
      "Hazard Ratio"    = HR,
      "CI95"            = CI,
      "P.value"         = Pvalue,
      stringsAsFactors  = FALSE
    )
    return(Unicox)
  }

  UniVar_list <- lapply(varNames, UniCox)
  UniVar <- ldply(UniVar_list, data.frame)

  multi_vars <- varNames
  if (filter_by_p_value) {
    sig_vars <- UniVar$Characteristics[!is.na(UniVar$P.value) &
                                         UniVar$P.value <= p_cutoff]
    if (length(sig_vars) == 0) {
      warning("没有单因素 p <", p_cutoff,
              " 的协变量，多因素 Cox 将不会拟合。")
      multi_vars <- character(0)
    } else {
      multi_vars <- sig_vars
    }
  }

  if (length(multi_vars) > 0) {
    fml <- as.formula(paste0("surv ~ ",
                             paste0(multi_vars, collapse = " + ")))
    multicox <- coxph(fml, data = data)
    multisum <- summary(multicox)
    
    multiname <- rownames(multisum$coefficients)  
    mCIL <- round(multisum$conf.int[, 3], 4)
    mCIU <- round(multisum$conf.int[, 4], 4)
    mCI <- paste0(mCIL, "-", mCIU)
    mPvalue <- round(multisum$coefficients[, 5], 4)
    mHR <- round(multisum$coefficients[, 2], 4)
    
    mulcox <- data.frame(
      "Characteristics"    = multiname,
      "Hazard Ratio.multi" = mHR,
      "CI95.multi"         = mCI,
      "P.value.multi"      = mPvalue,
      stringsAsFactors     = FALSE
    )
  } else {
    mulcox <- data.frame(
      "Characteristics"    = character(0),
      "Hazard Ratio.multi" = numeric(0),
      "CI95.multi"         = character(0),
      "P.value.multi"      = numeric(0),
      stringsAsFactors     = FALSE
    )
  }
  
  Final <- merge(UniVar, mulcox,
                 by = "Characteristics",
                 all = TRUE,
                 sort = FALSE) 
  Final$Characteristics <- factor(Final$Characteristics,
                                  levels = varNames)
  Final <- Final[order(Final$Characteristics), ]
  Final$Characteristics <- as.character(Final$Characteristics)
  
  return(Final)
}
epi_geo_meta <- read.csv('./Survival_analysis/Input/Cox_input/Epi/GEO-meta/merged.uni_multi_cox_input.csv', stringsAsFactors = F)
epi_gse26549 <- read.csv('./Survival_analysis/Input/Cox_input/Epi/GSE26549/merged.uni_multi_cox_input.csv', stringsAsFactors = F)
epi_tcga <- read.csv('./Survival_analysis/Input/Cox_input/Epi/TCGA/merged.uni_multi_cox_input.csv', stringsAsFactors = F)
out_1 <- cox_pipe(epi_geo_meta, var_n = 10, filter_by_p_value = F, p_cutoff = 0.05)
print(out_1)
write.csv(out_1, './Survival_analysis/Output/Epi_GEO-meta_COX_res.csv', row.names = FALSE, na = "")
out_2 <- cox_pipe(epi_gse26549, var_n = 11, filter_by_p_value = F, p_cutoff = 0.05)
print(out_2)
write.csv(out_2, './Survival_analysis/Output/Epi_GSE26549_COX_res.csv', row.names = FALSE, na = "")
out_3 <- cox_pipe(epi_tcga, var_n = 13, filter_by_p_value = F, p_cutoff = 0.05)
print(out_3)
write.csv(out_3, './Survival_analysis/Output/Epi_TCGA_COX_res.csv', row.names = FALSE, na = "")
fibro_geo_meta <- read.csv('./Survival_analysis/Input/Cox_input/Fibro/GEO-meta/merged.uni_multi_cox_input.csv', stringsAsFactors = F)
fibro_gse26549 <- read.csv('./Survival_analysis/Input/Cox_input/Fibro/GSE26549/merged.uni_multi_cox_input.csv', stringsAsFactors = F)
fibro_tcga <- read.csv('./Survival_analysis/Input/Cox_input/Fibro/TCGA/merged.uni_multi_cox_input.csv', stringsAsFactors = F)
fibro_out_1 <- cox_pipe(fibro_geo_meta, var_n = 9, filter_by_p_value = F, p_cutoff = 0.05)
print(fibro_out_1)
write.csv(fibro_out_1, './Survival_analysis/Output/Fibro_GEO-meta_COX_res.csv', row.names = FALSE, na = "")
fibro_out_2 <- cox_pipe(fibro_gse26549, var_n = 10, filter_by_p_value = F, p_cutoff = 0.05)
print(fibro_out_2)
write.csv(fibro_out_2, './Survival_analysis/Output/Fibro_GSE26549_COX_res.csv', row.names = FALSE, na = "")
fibro_out_3 <- cox_pipe(fibro_tcga, var_n = 12, filter_by_p_value = F, p_cutoff = 0.05)
print(fibro_out_3)
write.csv(fibro_out_3, './Survival_analysis/Output/Fibro_TCGA_COX_res.csv', row.names = FALSE, na = "")
plot_cox_forest <- function(cox_df,
                            out_dir = ".",
                            title_tag = "cox",
                            var_order = NULL  
                            ) {
  library(forestplot)
  library(stringr)
  library(grid)

  needed_cols <- c("Characteristics",
                   "Hazard.Ratio", "CI95", "P.value",
                   "Hazard.Ratio.multi", "CI95.multi", "P.value.multi")
  if (!all(needed_cols %in% colnames(cox_df))) {
    stop("cox_df 缺少必要列，请确认列名为：\n",
         paste(needed_cols, collapse = ", "))
  }

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  if (!is.null(var_order)) {
    not_in_df <- setdiff(var_order, cox_df$Characteristics)
    if (length(not_in_df) > 0) {
      warning("var_order 中这些变量不在 cox_df$Characteristics 中：",
              paste(not_in_df, collapse = ", "))
    }
    all_char <- c(var_order, setdiff(cox_df$Characteristics, var_order))
    cox_df <- cox_df[match(all_char, cox_df$Characteristics), ]
  } else {
    if ("Risk" %in% cox_df$Characteristics) {
      cox_df <- cox_df[order(cox_df$Characteristics == "Risk"), ]
    }
  }

  ci_uni   <- str_split_fixed(cox_df$CI95, "-", 2)
  ci_multi <- str_split_fixed(cox_df$CI95.multi, "-", 2)

  cox_df$Low.uni    <- as.numeric(ci_uni[, 1])
  cox_df$High.uni   <- as.numeric(ci_uni[, 2])
  cox_df$Low.multi  <- as.numeric(ci_multi[, 1])
  cox_df$High.multi <- as.numeric(ci_multi[, 2])

  .calc_xticks <- function(low, high) {
      vals <- c(low, high)
      vals <- vals[is.finite(vals)]
    
      rng <- range(vals, na.rm = TRUE)
      rng[1] <- max(floor(rng[1]), 0)
      rng[2] <- ceiling(rng[2])
    
      ticks <- seq(rng[1], rng[2], by = 1)
      if (!1 %in% ticks) ticks <- sort(c(ticks, 1))
      ticks
    }

  xticks_uni   <- .calc_xticks(cox_df$Low.uni,  cox_df$High.uni)
  xticks_multi <- .calc_xticks(cox_df$Low.multi, cox_df$High.multi)

  tabletext_uni <- cbind(
    c("Variable", cox_df$Characteristics),
    c(
      "HR (95% CI)",
      paste0(
        sprintf("%.2f", cox_df$Hazard.Ratio), " (",
        sprintf("%.2f", cox_df$Low.uni), "-",
        sprintf("%.2f", cox_df$High.uni), ")"
      )
    ),
    c(
      "P value",
      ifelse(cox_df$P.value < 0.001,
             "<0.001",
             sprintf("%.3f", cox_df$P.value))
    )
  )

  last_line_uni <- nrow(tabletext_uni) + 1
  hrzl_lines_uni <- list(
    "2" = gpar(lwd = 2, col = "black")
  )
  hrzl_lines_uni[[as.character(last_line_uni)]] <- gpar(lwd = 2, col = "black")

  tabletext_multi <- cbind(
    c("Variable", cox_df$Characteristics),
    c(
      "HR (95% CI)",
      paste0(
        sprintf("%.2f", cox_df$Hazard.Ratio.multi), " (",
        sprintf("%.2f", cox_df$Low.multi), "-",
        sprintf("%.2f", cox_df$High.multi), ")"
      )
    ),
    c(
      "P value",
      ifelse(cox_df$P.value.multi < 0.001,
             "<0.001",
             sprintf("%.3f", cox_df$P.value.multi))
    )
  )

  last_line_multi <- nrow(tabletext_multi) + 1
  hrzl_lines_multi <- list(
    "2" = gpar(lwd = 2, col = "black")
  )
  hrzl_lines_multi[[as.character(last_line_multi)]] <- gpar(lwd = 2, col = "black")

  pdf_file_uni <- file.path(out_dir, paste0("forest_uni_", title_tag, ".pdf"))
  pdf(pdf_file_uni, width = 12, height = 12,
      useDingbats = FALSE,
      family = "Helvetica")

  p_uni <- forestplot(
    labeltext = tabletext_uni,
    mean  = c(NA, cox_df$Hazard.Ratio),
    lower = c(NA, cox_df$Low.uni),
    upper = c(NA, cox_df$High.uni),

    graph.pos   = 2,
    graphwidth  = unit(.15, "npc"),
    fn.ci_norm  = "fpDrawCircleCI",
    col = fpColors(
      box  = "steelblue",
      lines = "black",
      zero  = "black"
    ),
    boxsize = 0.2,

    lwd.ci = 2,
    ci.vertices.height = 0.1,
    ci.vertices = TRUE,
    lwd.zero = 2,
    grid = structure(
      1, 
      gp = gpar(col = "black", lty = 2, lwd = 2)
    ),
    xticks   = xticks_uni,
    lwd.xaxis = 2,
    xlab = "Hazard Ratio (Univariate Cox)",

    hrzl_lines = hrzl_lines_uni,

    txt_gp = fpTxtGp(
      label = gpar(cex = 1),
      ticks = gpar(cex = 1),
      xlab  = gpar(cex = 1),
      title = gpar(cex = 1)
    ),
    lineheight = unit(0.75, "cm"),
    colgap     = unit(0, "cm"),
    mar        = unit(rep(1.25, times = 6), "cm"),
    new_page   = FALSE
  )
  print(p_uni)
  dev.off()

  pdf_file_multi <- file.path(out_dir, paste0("forest_multi_", title_tag, ".pdf"))
  pdf(pdf_file_multi, width = 12, height = 12,
      useDingbats = FALSE,
      family = "Helvetica")

  p_multi <- forestplot(
    labeltext = tabletext_multi,
    mean  = c(NA, cox_df$Hazard.Ratio.multi),
    lower = c(NA, cox_df$Low.multi),
    upper = c(NA, cox_df$High.multi),

    graph.pos   = 2,
    graphwidth  = unit(.15, "npc"),
    fn.ci_norm  = "fpDrawCircleCI",
    col = fpColors(
      box  = "steelblue",
      lines = "black",
      zero  = "black"
    ),
    boxsize = 0.2,

    lwd.ci = 2,
    ci.vertices.height = 0.1,
    ci.vertices = TRUE,
    lwd.zero = 2,
    grid = structure(
      1,
      gp = gpar(col = "black", lty = 2, lwd = 2)
    ),
    xticks   = xticks_multi,
    lwd.xaxis = 2,
    xlab = "Hazard Ratio (Multivariate Cox)",

    hrzl_lines = hrzl_lines_multi,

    txt_gp = fpTxtGp(
      label = gpar(cex = 1),
      ticks = gpar(cex = 1),
      xlab  = gpar(cex = 1),
      title = gpar(cex = 1)
    ),
    lineheight = unit(0.75, "cm"),
    colgap     = unit(0, "cm"),
    mar        = unit(rep(1.25, times = 6), "cm"),
    new_page   = FALSE
  )

  print(p_multi)
  dev.off()
}
print(out_1)
plot_cox_forest(out_1, out_dir = './Survival_analysis/Output/Epi/GEO-meta', 
                var_order = c('Stage_harmonized', 'Gender_harmonized', 'Age_harmonized', 'COL17A1_Epi', 'MKI67_Epi', 'GPX2_Epi', 'SBSN_Epi'))
print(out_2$Characteristics)
plot_cox_forest(out_2, out_dir = './Survival_analysis/Output/Epi/GSE26549', 
                var_order = c('Gender_harmonized', 'Age_harmonized', 'Smoking_harmonized', 'Alcohol_harmonized', 'COL17A1_Epi', 'MKI67_Epi', 'GPX2_Epi', 'SBSN_Epi'))
print(out_3$Characteristics)
plot_cox_forest(out_3, out_dir = './Survival_analysis/Output/Epi/TCGA', 
                var_order = c('Stage_harmonized', 'Gender_harmonized', 'Age_harmonized', 'Smoking_harmonized', 'Alcohol_harmonized', 'HPV_status_harmonized', 'COL17A1_Epi', 'MKI67_Epi', 'GPX2_Epi', 'SBSN_Epi'))
print(fibro_out_1$Characteristics)
plot_cox_forest(fibro_out_1, out_dir = './Survival_analysis/Output/Fibro/GEO-meta', 
                var_order = c('Stage_harmonized', 'Gender_harmonized', 'Age_harmonized', 'CFD_Fibro', 'CARBP1_Fibro', 'Myofibroblast'))
print(fibro_out_2$Characteristics)
plot_cox_forest(fibro_out_2, out_dir = './Survival_analysis/Output/Fibro/GSE26549', 
                var_order = c('Gender_harmonized', 'Age_harmonized', 'Smoking_harmonized', 'Alcohol_harmonized', 'CFD_Fibro', 'CARBP1_Fibro', 'Myofibroblast'))
print(fibro_out_3$Characteristics)
plot_cox_forest(fibro_out_3, out_dir = './Survival_analysis/Output/Fibro/TCGA', 
                var_order = c('Stage_harmonized', 'Gender_harmonized', 'Age_harmonized', 'Smoking_harmonized', 'Alcohol_harmonized', 'HPV_status_harmonized', 'CFD_Fibro', 'CARBP1_Fibro', 'Myofibroblast'))
