# Combine per-cohort DA results into one table per model.
#
# One row per cohort x target: the per-target linear-model estimate (logFC, SE,
# p) from Script 03 plus the QC annotations from Script 04. Written once per
# model so the meta-analysis can be run on each:
#   plate_adj  - per-cohort plate rule (covariate for ELPD, UMKLM_KUL (P118,
#                merged), TRAPCAF, CANDAS-SMPD; none for NTUH, where the
#                within-plate contrast would rest on a small subset of cases)
#   plate_none - no plate term in any cohort
#   sva        - phenotype-protected SVA in every cohort (written only for
#                cohorts whose _plate_sva folders exist)

library(tidyverse)

source("~/proteomics_nulisa/scripts/nulisa_pipeline_v2/quick_pass/00_config_v2.R")
BASE <- PROJECT_DIR

cohorts <- c(p111  = "ELPD",
             p118 = "UMKLM_KUL",
             p121  = "NTUH",
             p136  = "TRAPCAF",
             p143  = "Nigeria-PD",
             p149  = "CANDAS-SMPD")

plate_rule <- list(
  plate_adj  = c(p111 = "_plate_adj",  p118 = "_plate_adj",
  p121 = "_plate_none", p136  = "_plate_adj",  p143  = "_plate_none",
  p149 = "_plate_adj"),
  
plate_none = setNames(rep("_plate_none", length(cohorts)), names(cohorts)),
  sva        = setNames(rep("_plate_sva",  length(cohorts)), names(cohorts))
)

read_cohort <- function(co, label, sfx) {
  base <- file.path(BASE, paste0("results_", co, "_v2_filter_none"))
  da_f <- file.path(base, paste0("03_differential_abundance", sfx), "DA_results_all_Neuro220.csv")
  qc_f <- file.path(base, paste0("04_qc_report", sfx), "qc_target_Neuro220.csv")

  if (!file.exists(da_f)) { cat("  no DA results:", co, sfx, "\n"); return(NULL) }
  da <- read.csv(da_f, stringsAsFactors = FALSE)

  # QC adds plate ICC, CVs and the zero-imbalance flag.
  if (file.exists(qc_f)) {
    qc <- read.csv(qc_f, stringsAsFactors = FALSE)
    qc <- select(qc, Target, !any_of(names(da)))
    da <- left_join(da, qc, by = "Target")
    stopifnot(nrow(da) == nrow(read.csv(da_f)))   # join must not add rows
  } else {
    cat("  no QC file:", co, sfx, "- QC columns will be NA\n")
  }

  da %>% mutate(cohort = co, cohort_label = label, plate_model = sfx, .before = 1)
}

for (rule in names(plate_rule)) {
  sfx <- plate_rule[[rule]]
  combined <- bind_rows(lapply(names(cohorts), function(co)
    read_cohort(co, cohorts[[co]], sfx[[co]])))
  if (nrow(combined) == 0) next

  cat("\n", rule, "\n")
  print(count(combined, cohort_label, name = "n_targets"))
  write.csv(combined, file.path(BASE, paste0("DA_all_cohorts_combined_", rule, ".csv")),
            row.names = FALSE)
}