# Script 03 (v2): Differential Abundance (DA) - PD vs Control
# Description: DA for the cohort and panels set in 00_config_v2.R.
#              Annotates hits with target QC, vendor flags, and SC drift.




#  0. Setup 

library(SummarizedExperiment)
library(tidyverse)
library(limma)
library(ggrepel)

source("~/proteomics_nulisa/scripts/nulisa_pipeline_v2/quick_pass/00_config_v2.R")

dir.create(RESULTS_03_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(SC_DRIFT_FILE))
  stop("SC drift file not found - run 02b first:\n  ", SC_DRIFT_FILE)

# load META6 GWAS genes
meta6_genes  <- read.csv(META6_FILE, stringsAsFactors = FALSE)
meta6_unique <- unique(meta6_genes$Meta6_Genes)



#  1. Build the analysis cohort 

build_cohort <- function(se) {
  meta <- as.data.frame(colData(se))
  keep <- meta$phenotype_clean %in% c("PD", "Control")
  keep[is.na(keep)] <- FALSE
  se <- se[, keep]
  meta <- as.data.frame(colData(se))

  # exclude apoe sample swaps 
  n_before <- ncol(se)
  se <- se[, !colnames(se) %in% IDENTITY_EXCLUDE]
  if (ncol(se) < n_before) {
    cat("Excluded", n_before - ncol(se), "APOE sample-swap failure(s)\n")
    meta <- as.data.frame(colData(se))
  }

  # exclude severe low-detectability samples
  if (!is.null(SAMPLE_DETECT_EXCLUDE)) {
    bad <- !is.na(meta$pct_above_lod) & meta$pct_above_lod < SAMPLE_DETECT_EXCLUDE
    if (any(bad)) {
      cat("Excluding", sum(bad), "sample(s) below", SAMPLE_DETECT_EXCLUDE,
          "% detectability:\n")
      print(meta[bad, c("SampleName", "PlateID", "pct_above_lod")])
      se <- se[, !bad]
      meta <- as.data.frame(colData(se))
    }
  }

  # complete covariates
  cc <- complete.cases(meta[, c("phenotype_clean", DA_COVARIATES), drop = FALSE])
  se <- se[, cc]

  # targets need enough measured samples to model
  n_measured <- rowSums(!is.na(assay(se, "npq")))
  cat("Targets with fewer than", MIN_SAMPLES_PER_TARGET,
      "measured samples:", sum(n_measured < MIN_SAMPLES_PER_TARGET), "\n")
  se <- se[n_measured >= MIN_SAMPLES_PER_TARGET, ]

  colData(se)$phenotype_clean <- factor(colData(se)$phenotype_clean,
                                        levels = c("Control", "PD"))
  colData(se)$sex_clean <- factor(colData(se)$sex_clean)
  colData(se)$PlateID   <- droplevels(factor(colData(se)$PlateID))

  # APOE4 is a binary carrier readout so it is not modelled
  if (length(EXCLUDE_FROM_DA) > 0) {
    drop <- rownames(se) %in% EXCLUDE_FROM_DA
    if (any(drop)) {
      cat("Excluded from DA per config:", paste(rownames(se)[drop], collapse = ", "), "\n")
      se <- se[!drop, ]
    }
  }

  se
}



#  2. Fit the linear model per target

fit_da <- function(se) {
  meta <- as.data.frame(colData(se))
  expr <- assay(se, "npq")
  form <- paste("~ phenotype_clean +", paste(DA_COVARIATES, collapse = " + "))
  if (PLATE_CORRECTION == "covariate") form <- paste(form, "+ PlateID")
  design <- model.matrix(as.formula(form), data = meta)
  if (PLATE_CORRECTION == "sva") {
    if (!requireNamespace("sva", quietly = TRUE))
      stop("PLATE_CORRECTION is 'sva' but the sva package is not installed.")
    mod0 <- model.matrix(as.formula(paste("~", paste(DA_COVARIATES, collapse = " + "))),
                         data = meta)
    sv <- sva::sva(expr, design, mod0)
    if (sv$n.sv > 0) {
      design <- cbind(design, sv$sv)
      colnames(design)[(ncol(design) - sv$n.sv + 1):ncol(design)] <-
        paste0("SV", seq_len(sv$n.sv))
    }
    cat("  SVA estimated", sv$n.sv, "surrogate variable(s)\n")
  }
  stopifnot(nrow(design) == ncol(expr))   # a dropped row here means NA covariates
    
  fit <- lmFit(expr, design)

  coef_name <- "phenotype_cleanPD"

  se_ord <- fit$sigma * fit$stdev.unscaled[, coef_name]
  t_ord  <- fit$coefficients[, coef_name] / se_ord
  p_ord  <- 2 * pt(-abs(t_ord), df = fit$df.residual)

  tibble(Target    = rownames(expr),
         logFC     = fit$coefficients[, coef_name],
         SE        = se_ord,
         t         = t_ord,
         P.Value   = p_ord,
         adj.P.Val = p.adjust(p_ord, "BH"),
         df        = fit$df.residual,
         n         = ncol(expr)) %>%
    arrange(P.Value)
}



#  3. Genomic inflation  

lambda_from_p <- function(p) {
  p <- p[!is.na(p) & p > 0]
  median(qchisq(1 - p, 1)) / qchisq(0.5, 1)
}


#  4. PCA outlier detection

detect_outliers_pca <- function(expr_mat, meta_df, panel_name,
                                sd_threshold = OUTLIER_SD_THRESHOLD) {
  pca <- prcomp(t(expr_mat), scale. = TRUE, center = TRUE)
  var_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  
  # distance from centroid using first 5 PCs
  pc_scores <- pca$x[, 1:5]
  distances <- sqrt(rowSums(pc_scores^2))
  dist_threshold <- mean(distances) + sd_threshold * sd(distances)
  outlier_samples <- names(distances)[distances > dist_threshold]
  cat("Outliers found:", length(outlier_samples), "/", length(distances), "\n")
  # plot (outliers in black)
  pheno <- meta_df$phenotype_clean[match(names(distances), meta_df$SampleName)]
  point_colors <- ifelse(pheno == "PD", "steelblue", "tomato")
  point_colors[distances > dist_threshold] <- "black"
  png_file <- file.path(RESULTS_03_DIR, paste0("pca_outliers_", panel_name, ".png"))
  png(png_file, width = 10, height = 8, units = "in", res = 300)
  plot(pca$x[, 1], pca$x[, 2],
       col = point_colors, pch = 16, cex = 1.2,
       xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
       main = paste0("PCA Outlier Detection - ", panel_name),
       sub = paste0("Outliers (black) = >", sd_threshold, " SD beyond the mean distance from centroid (PCs 1-5)"))
  legend("topright", legend = c("PD", "Control", "Outlier"),
         col = c("steelblue", "tomato", "black"), pch = 16)
  dev.off()

  if (length(outlier_samples) > 0) {
    outlier_df <- data.frame(
      SampleName = outlier_samples,
      distance   = round(distances[outlier_samples], 3),
      threshold  = round(dist_threshold, 3)
    )
    write.csv(outlier_df,
              file.path(RESULTS_03_DIR, paste0("outliers_", panel_name, ".csv")),
              row.names = FALSE)
  } else {
    cat("No outliers detected for", panel_name, "\n")
  }
  return(outlier_samples)
}



# 5. Differential abundance for one panel 

run_differential_abundance <- function(se_path, panel_name) {
  cat("=====", panel_name, "=====\n")
  se_all <- readRDS(se_path)
  se <- build_cohort(se_all)
  meta <- as.data.frame(colData(se))
  cat("DA cohort:", ncol(se), "samples",
      "(PD:", sum(meta$phenotype_clean == "PD"),
      ", Control:", sum(meta$phenotype_clean == "Control"), ") |",
      nrow(se), "targets\n")
  stopifnot(!any(duplicated(meta$DONOR_ID)))
  
  # Plate balance. If plates are single-phenotype, plate-as-covariate removes
  # real signal rather than batch.
  plate_tab <- table(meta$PlateID, meta$phenotype_clean)
  n_single <- sum(rowSums(plate_tab > 0) < 2)
  cat("Plates:", nrow(plate_tab), "| single-phenotype plates:", n_single, "\n")
  if (n_single > 0) {
    cat("WARNING: some plates are single-phenotype;",
        "review before trusting the plate covariate.\n")
    print(plate_tab[rowSums(plate_tab > 0) < 2, , drop = FALSE])
  }
  outliers <- detect_outliers_pca(assay(se, "npq"), meta, panel_name)
  cat("  ", if (REMOVE_OUTLIERS) "REMOVED from the model"
            else "KEPT in the model; reported as a sensitivity column", "\n")
  if (REMOVE_OUTLIERS && length(outliers) > 0)
    se <- se[, !colnames(se) %in% outliers]


  #  main model
  results <- fit_da(se)
  lam <- lambda_from_p(results$P.Value)
  cat(sprintf("lambda = %.3f (reported only; not used to adjust p-values)\n", lam))
  
  #  sensitivity: drop SampleQC-warned samples 
  warned <- colnames(se)[colData(se)$sample_qc_flag %in% TRUE]
  sens_qc <- NULL
  if (length(warned) > 0 && length(warned) < ncol(se) - 20) {
    sens_qc <- fit_da(se[, !colnames(se) %in% warned]) %>%
      select(Target, logFC_noQCwarn = logFC, P_noQCwarn = P.Value)
    cat("Sensitivity: re-fit without", length(warned), "SampleQC-warned samples\n")
  }

  #  sensitivity: drop PCA outliers 
  sens_out <- NULL
  if (!REMOVE_OUTLIERS && length(outliers) > 0) {
    sens_out <- fit_da(se[, !colnames(se) %in% outliers]) %>%
      select(Target, logFC_noOutlier = logFC, P_noOutlier = P.Value)
    cat("Sensitivity: re-fit without", length(outliers), "PCA outliers\n")
  }

 # track non-zero data by phenotype 
   expr  <- assay(se, "npq")
  is_pd <- colData(se)$phenotype_clean == "PD"
  nonzero <- tibble(
    Target           = rownames(expr),
    pct_nonzero_pd   = round(100 * rowMeans(expr[, is_pd]  > 0), 1),
    pct_nonzero_ctrl = round(100 * rowMeans(expr[, !is_pd] > 0), 1))

  results <- results %>% left_join(nonzero, by = "Target")


  #  assemble 
  target_qc <- as.data.frame(rowData(se)) %>%
    select(Target, category, detectability, zero_prop, low_detect_flag)
  
  results <- results %>%
    mutate(
      panel      = panel_name,
      bonferroni = p.adjust(P.Value, method = "bonferroni"),
      in_meta6   = Target %in% meta6_unique,
      direction  = case_when(
        adj.P.Val < FDR_THRESHOLD & logFC > 0 ~ "Up in PD",
        adj.P.Val < FDR_THRESHOLD & logFC < 0 ~ "Down in PD",
        TRUE ~ "Not significant")
    ) %>%
    left_join(target_qc, by = "Target")
  if (!is.null(sens_qc))  results <- results %>% left_join(sens_qc,  by = "Target")
  if (!is.null(sens_out)) results <- results %>% left_join(sens_out, by = "Target")

  # vendor flags come from Alamar Neuro220 data sheet
  results <- results %>%
    mutate(
      vendor_high_cv      = Target %in% VENDOR_HIGH_CV_PLASMA[[panel_name]],
      vendor_low_detect   = Target %in% VENDOR_LOW_DETECT[[panel_name]],
      cross_reactive      = Target %in% names(CROSS_REACTIVE_TARGETS[[panel_name]]),
      xr_partner_on_panel = Target %in% XR_PARTNER_ON_PANEL[[panel_name]],
      hemolysis_sensitive = Target %in% names(VENDOR_HEMOLYSIS_SENSITIVE[[panel_name]])
    )

  # Below 50% detectability: still modelled, but flagged not reportable
  results <- results %>%
    mutate(reportable = category == "high_abundance" |
                        (!is.na(detectability) & detectability >= REPORT_MIN_DETECT))
  cat("Targets below the", REPORT_MIN_DETECT * 100, "% reporting gate:", sum(!results$reportable), "\n")

   # annotate with script 02b inter- and intra-plate CV 
  sc_cv <- read.csv(SC_DRIFT_FILE) %>%
  select(Target, sc_inter_cv = inter_cv, sc_max_drift = max_abs_drift)

  results <- results %>% left_join(sc_cv, by = "Target")

  cat("Targets with no SC CV (below LOD in pooled plasma):",
      sum(is.na(results$sc_inter_cv)), "\n")
  cat("Significant hits over Alamar's 15% SC CV benchmark:", 
      sum(results$adj.P.Val < FDR_THRESHOLD & results$sc_inter_cv > 15, 
      na.rm = TRUE), "\n")
  
  # Robust = significant in the main model and still nominally significant
  # in every sensitivity re-fit that ran.
  results$robust <- results$adj.P.Val < FDR_THRESHOLD
  if ("P_noQCwarn" %in% names(results))
    results$robust <- results$robust & results$P_noQCwarn < 0.05
  if ("P_noOutlier" %in% names(results))
    results$robust <- results$robust & results$P_noOutlier < 0.05
  
  cat("\nSignificant at FDR <", FDR_THRESHOLD, ":",
      sum(results$adj.P.Val < FDR_THRESHOLD),
      "| of those, reportable (detectability >=", REPORT_MIN_DETECT, "):",
      sum(results$adj.P.Val < FDR_THRESHOLD & results$reportable, na.rm = TRUE),
      "| robust to sensitivity re-fits:", sum(results$robust, na.rm = TRUE), "\n")

  cat("\nTop 15:\n")
  print(results %>% select(Target, logFC, P.Value, adj.P.Val, direction,
                 detectability, reportable, sc_inter_cv, vendor_high_cv, robust) %>% head(15))

  write.csv(results,
            file.path(RESULTS_03_DIR, paste0("DA_results_all_", panel_name, ".csv")),
            row.names = FALSE)

  write.csv(results %>% filter(adj.P.Val < FDR_THRESHOLD),
            file.path(RESULTS_03_DIR, paste0("DA_significant_", panel_name, ".csv")),
            row.names = FALSE)

  write.csv(results %>% filter(in_meta6) %>% arrange(adj.P.Val),
            file.path(RESULTS_03_DIR, paste0("DA_meta6_", panel_name, ".csv")),
            row.names = FALSE)

  attr(results, "lambda") <- lam
  attr(results, "n_samples") <- ncol(se)
  results
}



# 6. Volcano plot 

make_volcano_plot <- function(results, panel_name) {
  results <- results %>%
    mutate(direction_nominal = case_when(
      P.Value < 0.05 & logFC > 0 ~ "Up in PD",
      P.Value < 0.05 & logFC < 0 ~ "Down in PD",
      TRUE ~ "Not significant"))
  fdr_pvals <- results$P.Value[results$adj.P.Val < FDR_THRESHOLD]
  threshold_y_fdr <- if (length(fdr_pvals) > 0) -log10(max(fdr_pvals)) else NA
 top_targets <- results %>% filter(adj.P.Val < FDR_THRESHOLD) %>%
    arrange(P.Value) %>% head(15) %>% pull(Target)
  to_label <- results %>% filter(Target %in% top_targets | (in_meta6 & P.Value < 0.05))
  p <- ggplot(results, aes(logFC, -log10(P.Value), color = direction_nominal)) +
    geom_point(alpha = 0.6, size = 2.5) +
    # gray rings = <50% detectable
    geom_point(data = filter(results, low_detect_flag),
               shape = 1, size = 4, colour = "grey30", stroke = 0.6) +
    geom_point(data = filter(results, in_meta6),
               shape = 21, fill = "orange", color = "black", size = 3.5) +
    geom_text_repel(data = to_label, aes(label = Target), color = "black",
                    size = 3, max.overlaps = 30, box.padding = 0.5,
                    show.legend = FALSE) +
    scale_color_manual(values = c("Up in PD" = "tomato",
                                  "Down in PD" = "steelblue",
                                  "Not significant" = "gray80")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray40") +
    labs(title = paste0(panel_name, ": PD vs Control (", COHORT, ")"),
         subtitle = paste0("Orange = META6 | grey ring = <", DET_THRESH * 100,
                           "% detectable | dotted = p=0.05",
                           if (!is.na(threshold_y_fdr))
                             paste0(" | dashed = FDR<", FDR_THRESHOLD) else ""),
         x = "log2 Fold Change (PD vs Control)", y = "-log10(p-value)",
         color = "Direction (nominal p<0.05)") +
    theme_bw()
  if (!is.na(threshold_y_fdr))
    p <- p + geom_hline(yintercept = threshold_y_fdr, linetype = "dashed",
                        color = "gray20", linewidth = 0.7)
  ggsave(file.path(RESULTS_03_DIR, paste0("volcano_", panel_name, ".png")),
         p, width = 10, height = 8, units = "in")
  p
}



# 7. Run all panels 

SE_PATHS <- list(Neuro220 = SE_MERGED)
results <- list()
for (p in PANELS) {
  results[[p]] <- run_differential_abundance(SE_PATHS[[p]], p)
  make_volcano_plot(results[[p]], p)
}


# 8. Summary 

build_summary <- function(results, panel_name) {
  data.frame(
    panel           = panel_name,
    n_samples       = attr(results, "n_samples"),
    n_tested        = nrow(results),
    n_sig           = sum(results$adj.P.Val < FDR_THRESHOLD),
    n_up            = sum(results$direction == "Up in PD"),
    n_down          = sum(results$direction == "Down in PD"),
    n_sig_robust    = sum(results$robust, na.rm = TRUE),
    n_sig_lowdetect = sum(results$adj.P.Val < FDR_THRESHOLD &
                            results$low_detect_flag, na.rm = TRUE),
    n_meta6_total   = sum(results$in_meta6),
    n_meta6_sig     = sum(results$in_meta6 & results$adj.P.Val < FDR_THRESHOLD),
    n_bonf          = sum(results$bonferroni < 0.05),
    lambda          = round(attr(results, "lambda"), 3),
    n_sig_lowdetect_std  = sum(results$adj.P.Val < FDR_THRESHOLD &
                                 results$low_detect_flag &
                                 results$category == "standard", na.rm = TRUE),
    n_sig_lowdetect_rare = sum(results$adj.P.Val < FDR_THRESHOLD &
                                 results$low_detect_flag &
                                 results$category == "rare_case", na.rm = TRUE),
    n_sig_sc_cv_over15   = sum(results$adj.P.Val < FDR_THRESHOLD &
                                 results$sc_inter_cv > 15, na.rm = TRUE),
    n_sig_vendor_high_cv = sum(results$adj.P.Val < FDR_THRESHOLD &
                                 results$vendor_high_cv, na.rm = TRUE),
    n_sig_reportable = sum(results$adj.P.Val < FDR_THRESHOLD & results$reportable, na.rm = TRUE),
    plate_correction = PLATE_CORRECTION
  )
}

summary_df <- bind_rows(lapply(PANELS, function(p) build_summary(results[[p]], p)))

write.csv(summary_df,
          file.path(RESULTS_03_DIR, "differential_abundance_summary.csv"),
          row.names = FALSE)

print(summary_df)
