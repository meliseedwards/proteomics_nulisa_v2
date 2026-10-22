# Script 04 v2: Quality Control Tables + Figures
# Date: August 2026
# Description: Post-DA quality control for Neuro220 panel. 


# 0. Setup 

library(SummarizedExperiment)
library(tidyverse)
library(lme4)
library(readxl)

source("~/proteomics_nulisa/scripts/nulisa_pipeline_v2/quick_pass/00_config_v2.R")

dir.create(RESULTS_04_DIR, showWarnings = FALSE, recursive = TRUE)



# QC thresholds from the NULISAseqR user guide Table 2.1 
ICC_MAX        <- 10    
INTRA_CV_MAX   <- 10    
INTER_CV_MAX   <- 15    
RUN_DETECT_MIN <- 90   
BATCH_F_P      <- 0.01  
BATCH_TUKEY_P  <- 0.01  

# covariates screened against the PCs
CANDIDATE_COVARIATES <- c("phenotype_clean", "PlateID", "age", "sex_clean",
                          "pct_above_lod", "sample_qc_flag", "study_arm",
                          "age_at_diagnosis", "family_history_pd", "region_for_qc")

# function to save png
save_png <- function(plot, name, panel, w = 10, h = 7) {
  ggsave(file.path(RESULTS_04_DIR, paste0(name, "_", panel, ".png")),
         plot, width = w, height = h, dpi = 300)
}


# 1. Cohort 
# identical to build_cohort() in script 03.
prepare_cohort <- function(se) {
  meta <- as.data.frame(colData(se))
  keep <- meta$phenotype_clean %in% c("PD", "Control")
  keep[is.na(keep)] <- FALSE
  se <- se[, keep]
  meta <- as.data.frame(colData(se))

  # remove apoe sample swaps
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

  cc <- complete.cases(meta[, c("phenotype_clean", DA_COVARIATES), drop = FALSE])
  se <- se[, cc]
  n_measured <- rowSums(!is.na(assay(se, "npq")))
  se <- se[n_measured >= MIN_SAMPLES_PER_TARGET, ]

  if (length(EXCLUDE_FROM_DA) > 0)
    se <- se[!rownames(se) %in% EXCLUDE_FROM_DA, ]

  se
}


#  2. ICC + F-test + Tukey, per target
# The NulisaSeqR user guide defines a batch-effect target as ALL three of these:
#   run ICC > 10%, F-test unadjusted p < 0.01, Tukey-adjusted p < 0.01 for >=1 run
# Note: ICC is inflated when lack of plate randomization; exploratory here
icc_one_protein <- function(npq_values, plate) {
  d <- data.frame(npq = npq_values, plate = factor(plate))
  
  # mixed model: plate as random effect -> variance components for ICC
  mm <- lmer(npq ~ 1 + (1 | plate), data = d,
             control = lmerControl(calc.derivs = FALSE))
  vc <- as.data.frame(VarCorr(mm))
  plate_var <- vc$vcov[vc$grp == "plate"]
  resid_var <- vc$vcov[vc$grp == "Residual"]
  icc <- 100 * plate_var / (plate_var + resid_var)
  rm(mm, vc); # the mixed model is large - drop it before fitting the aov
  
  # F-test: plate as fixed effect in an ordinary lm
  aov_fit <- aov(npq ~ plate, data = d)
  f_p <- summary(aov_fit)[[1]][["Pr(>F)"]][1]
  
  # Tukey: smallest adjusted p across all plate pairs
  tk <- suppressWarnings(TukeyHSD(aov_fit)$plate)
  tukey_min_p <- if (is.null(tk)) NA_real_ else min(tk[, "p adj"], na.rm = TRUE)
  rm(aov_fit, tk)
  data.frame(ICC = icc, plate_F_pvalue = f_p, tukey_min_p = tukey_min_p)
}

# Distribution shape per target.
skewness <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(NA_real_)
  mean((x - mean(x))^3) / stats::sd(x)^3
}



# 3. CV from SC controls 
# Computed in 02b_SC_drift_v2.R: pooled replicates across plates, sd/mean on
# unlogged values, below-LOD excluded. Equivalent to NULISAseqR::interCV with
# useMean = FALSE, method = "count" (the settings Banner used in the TAP report.)
if (!file.exists(SC_DRIFT_FILE)) 
  stop("SC drift file not found - run 02b_SC_drift_v2.R first:\n  ", SC_DRIFT_FILE)

cv_tbl <- read.csv(SC_DRIFT_FILE) %>%
  select(Target, intra_cv = intra_cv_median, intra_cv_max, worst_plate, inter_cv)


# 4. Variance in a PC explained by a covariate 

pc_r2 <- function(pc, covariate) {
  d <- data.frame(pc = pc, x = covariate)
  d <- d[complete.cases(d), ]
  if (nrow(d) < 10 || length(unique(d$x)) < 2) return(NA_real_)
  tryCatch(summary(lm(pc ~ x, data = d))$r.squared, error = function(e) NA_real_)
}


# 5. Everything for one panel

run_qc_panel <- function(se_path, npq_file, panel_name, da_file) {

  cat("\n==========", panel_name, "==========\n")
  se   <- prepare_cohort(readRDS(se_path))
  expr <- assay(se, "npq")
  meta <- as.data.frame(colData(se))
  plate <- factor(meta$PlateID)
  cat("Cohort:", ncol(expr), "samples |", nrow(expr), "targets |",
      nlevels(plate), "plates\n")

  # long format, used by several blocks below
  long <- as.data.frame(expr) %>%
    rownames_to_column("Target") %>%
    pivot_longer(-Target, names_to = "SampleName", values_to = "NPQ") %>%
    left_join(meta %>% select(SampleName, PlateID, phenotype_clean),
              by = "SampleName")



  #  5a. ICC 

  cat("Computing ICC and Tukey tests...\n")
  icc_raw <- map_dfr(rownames(expr), function(tg) {
    res <- tryCatch(icc_one_protein(expr[tg, ], plate),
                    error = function(e) data.frame(ICC = NA_real_,
                                                   plate_F_pvalue = NA_real_,
                                                   tukey_min_p = NA_real_))
    res$Target <- tg
    res
  })

  #  5b. skewness per target 
  shape <- data.frame(
    Target    = rownames(expr),
    skewness  = round(apply(expr, 1, skewness), 3)
  )



  #  5c. Read the NPQ workbook
   npq_long <- bind_rows(lapply(file.path(NPQ_DIR, npq_file), function(f) {
    read_excel(f, sheet = 1, skip = NPQ_SKIP, na = "NA") %>%
      mutate(across(any_of(c("NPQ", "targetLOD_NPQ")), as.numeric))
  }))

   # Fix jumbled Aβ target names
  fix_mojibake <- function(v) {
    y <- iconv(v, from = "UTF-8", to = "latin1")
    ok <- !is.na(y)
    v[ok] <- y[ok]
    Encoding(v) <- "UTF-8"
    v
  }
  npq_long$Target <- fix_mojibake(npq_long$Target)

  # same aliases as Script 01
  if ("SampleMatrixType" %in% names(npq_long))
    npq_long <- dplyr::rename(npq_long, Biofluid = SampleMatrixType)
  if ("targetLOD_NPQ" %in% names(npq_long))
    npq_long <- dplyr::rename(npq_long, LOD = targetLOD_NPQ)

  npq_long <- npq_long %>% mutate(NPQ = as.numeric(NPQ), LOD = as.numeric(LOD))


  # 5d. DA results 
  da <- read.csv(da_file, stringsAsFactors = FALSE) %>%
    select(Target, logFC, P.Value, adj.P.Val, direction, in_meta6,
           any_of(c("P_noQCwarn", "P_noOutlier", "robust")))


  # 5e. qc_target.csv 
  qc_target <- as.data.frame(rowData(se)) %>%
    select(Target, category, detectability, zero_prop, na_prop, low_detect_flag) %>%
    left_join(icc_raw,   by = "Target") %>%
    left_join(shape,     by = "Target") %>%
    left_join(cv_tbl,    by = "Target") %>%
    left_join(da,        by = "Target") %>%
    mutate(
      #  individual flags. Each is one criterion, checkable on its own. 
      sig            = !is.na(adj.P.Val) & adj.P.Val < FDR_THRESHOLD,
      low_batch_icc  = !is.na(ICC) & ICC < ICC_MAX,
      good_intra_cv  = !is.na(intra_cv) & intra_cv < INTRA_CV_MAX,
      good_inter_cv  = !is.na(inter_cv) & inter_cv < INTER_CV_MAX,
      well_detected  = !is.na(detectability) & detectability >= DET_THRESH,
      cv_measurable  = !is.na(intra_cv) & !is.na(inter_cv),
      detect_measurable = !is.na(detectability),
      # user guide's 3-criterion batch-effect target definition
      batch_effect_target = !is.na(ICC) & ICC > ICC_MAX &
                            !is.na(plate_F_pvalue) & plate_F_pvalue < BATCH_F_P &
                            !is.na(tukey_min_p) & tukey_min_p < BATCH_TUKEY_P,
      #  composite. TRUE only if every criterion above is met. 
      high_confidence = sig & low_batch_icc & good_intra_cv & good_inter_cv &
                        well_detected) %>%
    arrange(desc(sig), P.Value)

  cat("Targets flagged as batch-effect (ICC>", ICC_MAX, "% & F p<", BATCH_F_P,
      " & Tukey p<", BATCH_TUKEY_P, "): ",
      sum(qc_target$batch_effect_target, na.rm = TRUE), " of ", nrow(qc_target),
      "\n", sep = "")

  cat("Significant DA targets:", sum(qc_target$sig),
      "| of those high_confidence:", sum(qc_target$high_confidence, na.rm = TRUE), "\n")

  if (any(qc_target$sig & !qc_target$high_confidence, na.rm = TRUE)) {
    cat("Significant but NOT high_confidence - check these before presenting/publishing:\n")
    print(qc_target %>%
            filter(sig, !high_confidence) %>%
                        select(Target, adj.P.Val, ICC, intra_cv, intra_cv_max, worst_plate,
                   inter_cv, detectability, low_batch_icc, skewness,
                   good_intra_cv, cv_measurable, good_inter_cv, well_detected))
  }

  write.csv(qc_target, file.path(RESULTS_04_DIR,
            paste0("qc_target_", panel_name, ".csv")), row.names = FALSE)



  # 5f. qc_target_by_plate.csv 
  target_by_plate <- npq_long %>%
    filter(SampleName %in% colnames(expr), Target %in% rownames(expr)) %>%
    select(SampleName, Target, PlateID, NPQ, LOD) %>%
    group_by(Target, PlateID) %>%
    summarise(n            = dplyr::n(),
              median_npq   = round(median(NPQ, na.rm = TRUE), 3),
              detectability = round(mean(NPQ > LOD, na.rm = TRUE), 3),
              n_zero       = sum(!is.na(NPQ) & NPQ == 0),
              .groups = "drop") %>%
    group_by(Target) %>%
    mutate(delta_median = round(median_npq - median(median_npq, na.rm = TRUE), 3)) %>%
    ungroup()
  rm(npq_long); invisible(gc(FALSE))   # ~380k rows, free it before the figures
  write.csv(target_by_plate, file.path(RESULTS_04_DIR,
            paste0("qc_target_by_plate_", panel_name, ".csv")), row.names = FALSE)


  #  5g. PCA 
  pca <- prcomp(t(expr), scale. = TRUE, center = TRUE)
  var_exp <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  n_pc <- min(10, ncol(pca$x))
  pc_df <- as.data.frame(pca$x[, 1:n_pc])
  pc_df$SampleName <- colnames(expr)
  pc_meta <- left_join(pc_df, meta, by = "SampleName")


  #  5h. qc_sample.csv 
  qc_sample <- meta %>%
    select(SampleName, DONOR_ID, PlateID, SampleQC, sample_qc_flag,
           pct_above_lod, detectability_flag, phenotype_clean,
           age, sex_clean) %>%
    left_join(pc_df %>% select(SampleName, PC1:PC5), by = "SampleName")
  write.csv(qc_sample, file.path(RESULTS_04_DIR,
            paste0("qc_sample_", panel_name, ".csv")), row.names = FALSE)


  #  5i. qc_plate.csv 
  # Includes the PD/Control counts: if phenotype is unbalanced across plates,
  # the plate ANOVA in qc_pc_variance is confounded and cannot be read as batch.
  run_detect <- target_by_plate %>%
    group_by(PlateID) %>%
    summarise(run_detectability = round(100 * mean(detectability > 0.5, na.rm = TRUE), 1),
              .groups = "drop")

  qc_plate <- meta %>%
    group_by(PlateID) %>%
    summarise(n_samples   = dplyr::n(),
              n_PD        = sum(phenotype_clean == "PD"),
              n_Control   = sum(phenotype_clean == "Control"),
              pct_PD      = round(100 * mean(phenotype_clean == "PD"), 1),
              n_qc_warn   = sum(sample_qc_flag, na.rm = TRUE),
              median_sample_detect = round(median(pct_above_lod, na.rm = TRUE), 1),
              median_age  = round(median(age, na.rm = TRUE), 1),
              .groups = "drop") %>%
    left_join(run_detect, by = "PlateID") %>%
    mutate(run_detect_pass = run_detectability >= RUN_DETECT_MIN)

  # is phenotype balanced across plates
  bal_tab  <- table(meta$PlateID, meta$phenotype_clean)
  bal_test <- suppressWarnings(chisq.test(bal_tab))
  cat("Phenotype x plate balance: chi-sq p =", signif(bal_test$p.value, 3),
      if (bal_test$p.value < 0.05)
        "\n  -> UNBALANCED. The plate ANOVA below partly reflects phenotype, not batch\n"
      else "\n  -> balanced; plate effects are interpretable as batch.\n")
  qc_plate$phenotype_plate_chisq_p <- signif(bal_test$p.value, 4)
  cat("Run detectability (% of targets detectable on a plate): median ",
      round(median(qc_plate$run_detectability, na.rm = TRUE), 1),
      "%, plates below the guide's ", RUN_DETECT_MIN, "%: ",
      sum(!qc_plate$run_detect_pass), " of ", nrow(qc_plate), "\n", sep = "")
  if (sum(!qc_plate$run_detect_pass) > nrow(qc_plate) / 2)
    cat("Note: kept low-detectability targets lower this metric, not a plate failure\n")
  write.csv(qc_plate, file.path(RESULTS_04_DIR,
            paste0("qc_plate_", panel_name, ".csv")), row.names = FALSE)


  #  5j. covariate screen 
  # R-squared of each candidate covariate against each PC
  covars <- intersect(CANDIDATE_COVARIATES, colnames(pc_meta))
  cov_r2 <- expand.grid(covariate = covars, PC = paste0("PC", 1:n_pc),
                        stringsAsFactors = FALSE) %>%
    mutate(r2 = round(map2_dbl(PC, covariate,
                    ~ pc_r2(pc_meta[[.x]], pc_meta[[.y]])), 4))
  write.csv(cov_r2, file.path(RESULTS_04_DIR,
            paste0("qc_covariate_pc_r2_", panel_name, ".csv")), row.names = FALSE)

  top_cov <- cov_r2 %>%
    filter(PC %in% c("PC1", "PC2", "PC3")) %>%
    group_by(covariate) %>%
    summarise(max_r2 = ifelse(all(is.na(r2)), NA_real_, max(r2, na.rm = TRUE)),
              .groups = "drop") %>%
    arrange(desc(max_r2))
  cat("\nVariance in PC1-3 explained by each candidate covariate:\n")
  print(as.data.frame(top_cov))
  not_modelled <- top_cov %>%
    filter(!is.na(max_r2), max_r2 > 0.05, !covariate %in% c(DA_COVARIATES, "phenotype_clean",
                                            "PlateID"))

  if (nrow(not_modelled) > 0) {
    cat("Explains >5% of PC1-3 variance and is NOT in the model - consider:\n")
    print(as.data.frame(not_modelled))
    if ("pct_above_lod" %in% not_modelled$covariate)
      cat("  NOTE: pct_above_lod is calculated from this NPQ matrix, so it will\n",
          "  always track PC1 (overall signal level). It is a QC metric.\n")
  }


  #  5k. Zero report 
  ex  <- assay(se, "npq")
  qc <- qc_target
  pd  <- meta$phenotype_clean == "PD"
  ct  <- meta$phenotype_clean == "Control"

zero_report <- lapply(rownames(ex), function(t) {
  y <- ex[t, ]
  data.frame(
    Target        = t,
    n             = sum(!is.na(y)),
    pct_zero      = round(100 * mean(y == 0, na.rm = TRUE), 2),
    pct_zero_PD   = round(100 * mean(y[pd] == 0, na.rm = TRUE), 2),
    pct_zero_Ctrl = round(100 * mean(y[ct] == 0, na.rm = TRUE), 2))
}) %>% bind_rows() %>%
  mutate(zero_gap = round(abs(pct_zero_PD - pct_zero_Ctrl), 2)) %>%
  left_join(qc %>% select(Target, category, detectability, sig, adj.P.Val,
                          logFC, high_confidence, ICC, intra_cv, inter_cv),
            by = "Target") %>%
  mutate(zero_flag = pct_zero > 5,
         zero_imbalance_flag = zero_gap > 5,
         cohort = COHORT, panel = panel_name) %>%
  select(cohort, panel, Target, category, n, pct_zero, pct_zero_PD, pct_zero_Ctrl,
         zero_gap, zero_flag, zero_imbalance_flag, detectability,
         sig, logFC, adj.P.Val, high_confidence, ICC, intra_cv, inter_cv) %>%
  arrange(desc(zero_gap), desc(pct_zero))

write.csv(zero_report,
          file.path(RESULTS_04_DIR, paste0("zero_report_", panel_name, ".csv")),
          row.names = FALSE)

# put the zero columns on qc_target too, so one file has every flag
qc_target <- qc_target %>%
  left_join(zero_report %>% select(Target, zero_gap, zero_imbalance_flag),
            by = "Target")

write.csv(qc_target, file.path(RESULTS_04_DIR,
          paste0("qc_target_", panel_name, ".csv")), row.names = FALSE)

cat("targets with any exact zeros:", sum(zero_report$pct_zero > 0), "of", nrow(zero_report), "\n")
cat("targets >5% zeros:", sum(zero_report$zero_flag), "\n")
cat("targets with >5 point PD/Control zero gap:", sum(zero_report$zero_imbalance_flag), "\n\n")
print(head(zero_report %>% filter(pct_zero > 0) %>%
             select(Target, pct_zero, pct_zero_PD, pct_zero_Ctrl, zero_gap,
                    detectability, sig), 25))



  # =========================== FIGURES =====================================

  theme_set(theme_bw(base_size = 11))

  # 01. target x plate heatmap 
  hm_df <- target_by_plate %>%
    left_join(qc_target %>% select(Target, ICC), by = "Target") %>%
    mutate(Target = reorder(Target, ICC))

  p01 <- ggplot(hm_df, aes(PlateID, Target, fill = delta_median)) +
    geom_tile() +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato",
                         midpoint = 0, name = "median NPQ\n- target median") +
    labs(title = paste0(panel_name, ": target shift by plate"),
         subtitle = "Targets ordered by ICC (highest at top). A coloured stripe = plate effect.",
         x = "Plate", y = NULL) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
          axis.text.y = element_text(size = 4))
  save_png(p01, "01_target_by_plate_heatmap_all", panel_name, 11, 12)


  top_icc <- qc_target %>% arrange(desc(ICC)) %>% head(30) %>% pull(Target)

  p01b <- ggplot(hm_df %>% filter(Target %in% top_icc),
                 aes(PlateID, Target, fill = delta_median)) +
    geom_tile() +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato",
                         midpoint = 0, name = "median NPQ\n- target median") +
    labs(title = paste0(panel_name, ": 30 targets with the largest plate effects"),
         x = "Plate", y = NULL) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))

  save_png(p01b, "01b_target_by_plate_heatmap_top30", panel_name, 10, 8)



  # 02 / 03. Q-Q plots 
  # Observed NPQ against normal theoretical quantiles, coloured by plate, with
  # a reference line through the quartiles.
  make_qq <- function(targets, title, subtitle, fname) {
    d <- long %>%
      filter(Target %in% targets, !is.na(NPQ)) %>%
      group_by(Target) %>%
      mutate(theoretical = qnorm(ppoints(dplyr::n()))[rank(NPQ, ties.method = "first")]) %>%
      ungroup()
    refs <- d %>%
      group_by(Target) %>%
      summarise(y1 = quantile(NPQ, .25), y2 = quantile(NPQ, .75),
                x1 = qnorm(.25), x2 = qnorm(.75), .groups = "drop") %>%
      mutate(slope = (y2 - y1) / (x2 - x1), intercept = y1 - slope * x1)
    p <- ggplot(d, aes(theoretical, NPQ)) +
      geom_point(aes(colour = PlateID), size = 0.6, alpha = 0.7) +
      geom_abline(data = refs, aes(slope = slope, intercept = intercept),
                  colour = "red", linewidth = 0.5) +
      facet_wrap(~ Target, scales = "free_y", ncol = 3) +
      labs(title = title, subtitle = subtitle,
           x = "Theoretical quantiles", y = "NPQ") +
      theme(legend.position = "none")
    save_png(p, fname, panel_name, 10, 9)
  }
  hits <- qc_target %>% filter(sig) %>% arrange(P.Value) %>% head(9) %>% pull(Target)
  if (length(hits) > 0)
    make_qq(hits, paste0(panel_name, ": Q-Q for the top DA hits"),
            "Points coloured by plate; red line through the quartiles",
            "02_qq_top_hits")
  worst <- qc_target %>% arrange(detectability) %>% head(9) %>% pull(Target)
  make_qq(worst, paste0(panel_name, ": Q-Q for the 9 least detectable targets"),
          "A spike at the floor means the limma t-statistic is not measuring what it assumes",
          "03_qq_worst_detectability")



  # 04. PCA, 3 panels 
   pc_long <- pc_meta %>%
    select(PC1, PC2, Plate = PlateID, Phenotype = phenotype_clean,
           Sex = sex_clean) %>%
    pivot_longer(c(Plate, Phenotype, Sex), names_to = "variable", values_to = "value")

  p04 <- ggplot(pc_long, aes(PC1, PC2, colour = value)) +
    geom_point(size = 1.6, alpha = 0.7) +
    facet_wrap(~ variable) +
    labs(x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
         y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
         title = paste0(panel_name, ": PCA"), colour = NULL)

  save_png(p04, "04_pca_panels", panel_name, 12, 5)



  # 05. sample boxplots by plate
  p05 <- ggplot(long, aes(SampleName, NPQ, fill = PlateID)) +
    geom_boxplot(outlier.size = 0.2, linewidth = 0.2) +
    facet_grid(~ PlateID, scales = "free_x", space = "free_x") +
    labs(title = paste0(panel_name, ": NPQ distribution per sample, grouped by plate"),
         subtitle = "After normalisation these should sit at similar medians",
         x = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none",
          strip.text.x = element_text(size = 4, angle = 90))
  save_png(p05, "05_sample_boxplots_by_plate", panel_name, 16, 6)




# 06. sample detectability by plate 

  p06 <- ggplot(meta, aes(PlateID, pct_above_lod)) +
    geom_boxplot(outlier.size = 0.4, fill = "grey90") +
    geom_hline(yintercept = SAMPLE_DET_MIN * 100, linetype = "dashed",
               colour = "tomato") +
    labs(title = paste0(panel_name, ": sample detectability by plate"),
         subtitle = paste0("Dashed line = Alamar ", MATRIX_TYPE, " minimum (",
                           SAMPLE_DET_MIN * 100, "%)"),
         x = "Plate", y = "% of targets above LOD") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
  save_png(p06, "06_sample_detect_by_plate", panel_name, 11, 6)



  # 07. ICC distribution 
  # Show significant targets ordered by ICC.
  icc_lab <- qc_target %>% filter(sig) %>% arrange(desc(ICC)) %>% head(20)
  fig_title <- paste0(panel_name, ": DA hits ordered by plate ICC")
  if (nrow(icc_lab) == 0) {
    icc_lab <- qc_target %>% arrange(desc(ICC)) %>% head(12)
    fig_title <- paste0(panel_name, ": 12 most batch-sensitive targets (no DA hits)")
  }

  p07 <- ggplot(qc_target, aes(ICC)) +
    geom_histogram(binwidth = 1, fill = "grey70", colour = "white") +
    geom_vline(xintercept = ICC_MAX, linetype = "dashed", colour = "tomato") +
    labs(title = paste0(panel_name, ": plate ICC per target"),
         subtitle = paste0("Dashed line = user guide threshold (", ICC_MAX,
                           "%). Above it = substantial batch effect."),
         x = "ICC (% of variance from plate)", y = "Targets")
  save_png(p07, "07_icc_distribution", panel_name, 9, 5)

  p07b <- ggplot(icc_lab, aes(reorder(Target, ICC),
    ICC, fill = ICC > ICC_MAX)) +
    geom_col() + coord_flip() +
    geom_hline(yintercept = ICC_MAX, linetype = "dashed", colour = "tomato") +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "tomato"),
                      name = paste0("ICC > ", ICC_MAX, "%")) +
    labs(title = fig_title,
         subtitle = "Red = plate effect; treat hits with caution",
         x = NULL, y = "ICC (%)")
  save_png(p07b, "07b_icc_top_targets", panel_name, 8, 5)



  # 08. CV vs detectability 
  cv_long <- qc_target %>%
    select(Target, detectability, intra_cv, inter_cv, sig) %>%
    pivot_longer(c(intra_cv, inter_cv), names_to = "cv_type", values_to = "cv")
  p08 <- ggplot(cv_long, aes(detectability * 100, cv, colour = sig)) +
    geom_point(size = 1.4, alpha = 0.7) +
    facet_wrap(~ cv_type, scales = "free_y") +
    scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "tomato"),
                        name = "DA significant") +
    labs(title = paste0(panel_name, ": CV vs detectability"),
         subtitle = "CV is measured on SC wells, which are pooled PLASMA controls - not samples",
         x = "Target detectability (%)", y = "CV (%)")
  save_png(p08, "08_cv_vs_detectability", panel_name, 10, 5)


  # 09. covariate x PC heatmap 
  p09 <- ggplot(cov_r2, aes(factor(PC, levels = paste0("PC", 1:n_pc)),
                            reorder(covariate, r2, function(z) max(z, na.rm = TRUE)),
                            fill = r2)) +
    geom_tile() +
    geom_text(aes(label = ifelse(r2 > 0.05, sprintf("%.2f", r2), "")), size = 2.6) +
    scale_fill_gradient(low = "white", high = "darkred", name = "R-squared",
                        limits = c(0, 1)) +
    labs(title = paste0(panel_name, ": variance in each PC explained by each covariate"),
         subtitle = "Anything high here that is not already in the DA model is a candidate covariate",
         x = NULL, y = NULL)
  save_png(p09, "09_covariate_pc_r2", panel_name, 10, 6)


  # 10. sample correlation heatmap 
  set.seed(1)
  sub_n  <- min(150, ncol(expr))
  sub_id <- sample(colnames(expr), sub_n)
  cm <- cor(expr[, sub_id], use = "pairwise.complete.obs")
  ord <- hclust(as.dist(1 - cm))$order
  cm_df <- as.data.frame(cm[ord, ord]) %>%
    rownames_to_column("s1") %>%
    pivot_longer(-s1, names_to = "s2", values_to = "r") %>%
    mutate(s1 = factor(s1, levels = rownames(cm)[ord]),
           s2 = factor(s2, levels = rownames(cm)[ord]))
  p11 <- ggplot(cm_df, aes(s1, s2, fill = r)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Pearson r") +
    labs(title = paste0(panel_name, ": sample-sample correlation (", sub_n,
                        " random samples, clustered)"),
         subtitle = "Blocks on the diagonal = groups of similar samples",
         x = NULL, y = NULL) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  save_png(p11, "10_sample_correlation", panel_name, 8, 7)


  # 11. violins for significant targets passing all QC criteria 
  hc <- qc_target %>% filter(sig, high_confidence) %>% arrange(adj.P.Val) %>% head(12)
  if (nrow(hc) > 0) {
    vd <- as.data.frame(expr[hc$Target, , drop = FALSE]) %>%
      rownames_to_column("Target") %>%
      pivot_longer(-Target, names_to = "SampleName", values_to = "NPQ") %>%
      left_join(meta %>% select(SampleName, phenotype_clean), by = "SampleName") %>%
      left_join(hc %>% select(Target, logFC, adj.P.Val), by = "Target") %>%
      mutate(lab = paste0(Target, "  logFC ", sprintf("%+.2f", logFC),
                          ", FDR ", signif(adj.P.Val, 2)),
             phenotype_clean = factor(phenotype_clean, levels = c("Control", "PD")))

    p12 <- ggplot(vd, aes(phenotype_clean, NPQ, fill = phenotype_clean)) +
      geom_violin(alpha = 0.45, trim = FALSE) +
      geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white") +
      geom_jitter(width = 0.09, alpha = 0.3, size = 0.9) +
      facet_wrap(~ lab, scales = "free_y") +
      scale_fill_manual(values = c(Control = "tomato", PD = "steelblue")) +
      labs(title = paste0(toupper(COHORT), ": significant targets passing all QC criteria"),
           subtitle = paste0(sum(meta$phenotype_clean == "PD"), " PD vs ",
                             sum(meta$phenotype_clean == "Control"), " Control | plate ",
                             PLATE_CORRECTION),
           x = NULL, y = "NPQ (log2)") +
      theme(legend.position = "none")
    save_png(p12, "11_violin_high_confidence", panel_name, 11, 8)
  } else {
    cat("No high-confidence significant targets - violin figure skipped\n")
  }

  cat("Done:", panel_name, "\n")
  invisible(qc_target)
}



# 6. Run all panels 
SE_PATHS  <- list(Neuro220 = SE_MERGED)
NPQ_FILES <- list(Neuro220 = NPQ_FILE)

qc <- list()
for (p in PANELS) {
  qc[[p]] <- run_qc_panel(
    se_path  = SE_PATHS[[p]],
    npq_file   = NPQ_FILES[[p]],
    panel_name = p,
    da_file    = file.path(RESULTS_03_DIR, paste0("DA_results_all_", p, ".csv")))
}

# 7. Summary 
qc_summary <- bind_rows(lapply(PANELS, function(p) data.frame(
  panel             = p,
  n_targets         = nrow(qc[[p]]),
  n_sig             = sum(qc[[p]]$sig, na.rm = TRUE),
  n_high_confidence = sum(qc[[p]]$high_confidence, na.rm = TRUE),
  n_batch_effect    = sum(qc[[p]]$batch_effect_target, na.rm = TRUE),
  n_icc_over_thr    = sum(qc[[p]]$ICC > ICC_MAX, na.rm = TRUE),
  n_intra_cv_fail   = sum(qc[[p]]$intra_cv >= INTRA_CV_MAX, na.rm = TRUE),
  n_inter_cv_fail   = sum(qc[[p]]$inter_cv >= INTER_CV_MAX, na.rm = TRUE),
  n_cv_unmeasurable = sum(!qc[[p]]$cv_measurable),
  n_low_detect        = sum(qc[[p]]$detectability < DET_THRESH, na.rm = TRUE),
  n_zero_imbalance  = sum(qc[[p]]$zero_imbalance_flag, na.rm = TRUE),
  n_detect_unmeasurable = sum(!qc[[p]]$detect_measurable))))

write.csv(qc_summary, file.path(RESULTS_04_DIR, "qc_summary.csv"), row.names = FALSE)
print(qc_summary)