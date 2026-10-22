# Script 02 (v2): Merge GP2 Master Key (R12) with NULISA Data
# Date: October 2026
# Description: Loads SE objects from Script 01, joins the GP2 master key
#              using cohort donor ID, adds clinical metadata (diagnosis,
#              age, sex, ancestry) to colData, harmonizes phenotype to
#              PD/Control, and saves updated SE objects.


# 0. Setup
library(SummarizedExperiment)
library(tidyverse)
library(readxl)

source("~/proteomics_nulisa/scripts/nulisa_pipeline_v2/quick_pass/00_config_v2.R")

dir.create(RESULTS_02_DIR, showWarnings = FALSE, recursive = TRUE)


# 1. Load master key 
master_key <- read.csv(file.path(DATA_DIR, MASTER_KEY_FILE), stringsAsFactors = FALSE)

# Filter to study only
study_key <- master_key %>%
  filter(study %in% STUDY_FILTER) %>%
  mutate(clinical_id = as.character(clinical_id))
cat("Study samples in master key:", nrow(study_key), "\n")
stopifnot(nrow(study_key) > 0)

# Select only columns we need for downstream analysis
study_key_clean <- study_key %>%
  select(clinical_id, GP2ID, GP2sampleID, study,
         diagnosis, GP2_PHENO, GP2_phenotype, study_arm, study_type,
         biological_sex_for_qc, age, age_of_onset, age_at_diagnosis,
         age_at_death, age_at_last_follow_up, GP2_phenotype_for_qc,
         family_history_pd, nba_label, nba, wgs_label,
         race_for_qc, region_for_qc) %>%
mutate(join_key = .data[[GP2_JOIN_KEY]]) %>%
  filter(!is.na(join_key), trimws(join_key) != "")
cat("Master-key rows with a usable join key:", nrow(study_key_clean), "\n")

# Check for duplicate IDs
duplicated_ids <- study_key_clean %>%
  dplyr::count(join_key) %>%
  dplyr::filter(n > 1)
cat("Duplicate join keys in join key:", nrow(duplicated_ids), "\n")

if (nrow(duplicated_ids) > 0) {
  print(head(duplicated_ids))
  stop("Duplicate id(s) in the master key - the join would add rows.")
}


# 2. Merge metadata into SE
merge_metadata <- function(se_path, panel_name) {

  # load SE object from Script 01
  se <- readRDS(se_path)

  # extract current colData; DONOR_ID is the join key
  current_meta <- as.data.frame(colData(se)) %>%
    mutate(DONOR_ID = as.character(DONOR_ID))

  # join with master key 
   merged_meta <- current_meta %>%
    left_join(study_key_clean, by = c("DONOR_ID" = "join_key"))

  # Manifest fallback: some donors are not yet in R12, but the manifest carries
  # the same phenotype/age/sex fields. R12 wins where both exist.
  if (!is.null(MANIFEST_FILES)) {
      man <- bind_rows(lapply(seq_along(MANIFEST_FILES), function(i) {
      m   <- read.csv(file.path(MANIFEST_DIR, MANIFEST_FILES[i]), stringsAsFactors = FALSE)
      col <- if (length(MANIFEST_SAMPLE_COL) > 1) MANIFEST_SAMPLE_COL[i] else MANIFEST_SAMPLE_COL
      m %>% transmute(SampleName    = trimws(gsub("'", "", as.character(.data[[col]]))),
                      man_phenotype = GP2_phenotype,
                      man_age       = age,
                      man_sex       = biological_sex_for_qc)
    })) %>%
      distinct(SampleName, .keep_all = TRUE)

    merged_meta <- merged_meta %>%
      left_join(man, by = "SampleName") %>%
      mutate(across(c(GP2_phenotype, biological_sex_for_qc),
                    ~ ifelse(trimws(.x) == "", NA_character_, .x)),
             metadata_source = ifelse(is.na(GP2ID), "manifest", "master_key"),
             GP2_phenotype = coalesce(GP2_phenotype, man_phenotype),
             age           = coalesce(age, man_age),
             biological_sex_for_qc = coalesce(biological_sex_for_qc, man_sex))

    cat("Metadata source:\n"); print(table(merged_meta$metadata_source))
  }

  stopifnot(nrow(merged_meta) == ncol(se))

  # Report match success
  matched   <- sum(!is.na(merged_meta$GP2ID))
  unmatched <- nrow(merged_meta) - matched
  match_pct <- round(100 * matched / nrow(merged_meta), 1)

  cat("Samples matched to GP2:", matched, "/", nrow(merged_meta),
      "(", match_pct, "%)\n", sep = "")
  cat("Unmatched samples:", unmatched, "\n")
  cat("\nDiagnosis distribution (matched samples):\n")
  print(table(merged_meta$diagnosis, useNA = "always"))
  cat("\nGP2_PHENO distribution:\n")
  print(table(merged_meta$GP2_PHENO, useNA = "always"))
  if (sum(!is.na(merged_meta$age)) > 0) {
    cat("\nAge distribution:\n")
    print(summary(merged_meta$age))
  }


  #  Harmonize phenotype and sex labels 
  merged_meta <- merged_meta %>%
    mutate(
      sex_clean = case_when(
       trimws(biological_sex_for_qc) %in% c("F", "Female", "1") ~ "Female",
       trimws(biological_sex_for_qc) %in% c("M", "Male", "2")   ~ "Male",
       TRUE ~ NA_character_),
      # PD vs Control from standardized GP2 fields (prodromal/other dropped)
      phenotype_clean = case_when(
       trimws(GP2_phenotype) == "PD"      ~ "PD",
       trimws(GP2_phenotype) == "Control" ~ "Control",
       TRUE ~ NA_character_),
      # Unified ancestry: prefer NBA label; fall back to WGS label; else NA.
      ancestry = case_when(
      !is.na(nba_label) & nba_label != "" ~ nba_label,
      !is.na(wgs_label) & wgs_label != "" ~ wgs_label,
      TRUE                                ~ NA_character_
    ),
      ancestry_source = case_when(
      !is.na(nba_label) & nba_label != "" ~ "NBA",
      !is.na(wgs_label) & wgs_label != "" ~ "WGS",
      TRUE                                ~ NA_character_
    )
  )

  # exclude prodromal
  prodromal <- !is.na(merged_meta$diagnosis) & merged_meta$diagnosis %in% EXCLUDE_DIAGNOSES
  merged_meta$phenotype_clean[prodromal] <- NA

  # matched or manifest-sourced but GP2_phenotype was not PD or Control
  excluded <- merged_meta %>% filter(is.na(phenotype_clean))
  cat("Excluded (no clean PD/Control phenotype):", nrow(excluded), "\n")
  if (nrow(excluded) > 0) {
    print(excluded %>% dplyr::count(across(any_of(c("metadata_source", "GP2_phenotype", "diagnosis"))), sort = TRUE))
  }

  cat("\nHarmonized phenotype distribution:\n")
  print(table(merged_meta$phenotype_clean, useNA = "always"))

  cat("\nAncestry (unified: NBA then WGS) distribution:\n")
  print(table(merged_meta$ancestry, useNA = "always"))
  cat("\nAncestry source (NBA vs WGS):\n")
  print(table(merged_meta$ancestry_source, useNA = "always"))

  # Critical: merged_meta MUST be in same order as colnames(se)
  stopifnot(identical(as.character(merged_meta$SampleName), colnames(se)))
  rownames(merged_meta) <- merged_meta$SampleName
  colData(se) <- DataFrame(merged_meta)

# save updated SE
out_path <- file.path(RESULTS_02_DIR,
                        paste0("se_", COHORT, "_", panel_name, "_with_metadata.rds"))
saveRDS(se, out_path)
  cat("Saved:", basename(out_path), "\n")


# Donors not in R12. For manifest cohorts their phenotype/age/sex came from
  # the manifest instead, so they are still analysed - this list is for
  # tracking which donors still need a genotyping release (ancestry, variants).
  unmatched_ids <- merged_meta %>%
    filter(is.na(GP2ID)) %>%
    select(SampleName, DONOR_ID, PlateID, SampleQC, pct_above_lod)
  if (nrow(unmatched_ids) > 0) {
    write.csv(unmatched_ids,
              file.path(RESULTS_02_DIR,
                        paste0("not_in_r12_", COHORT, "_", panel_name, ".csv")),
              row.names = FALSE)
    cat("Wrote ", nrow(unmatched_ids),
        " donors not in R12 to not_in_r12_", COHORT, "_", panel_name, ".csv",
        " (metadata from manifest where available)\n", sep = "")
  }


  # Return summary stats
  return(list(
    panel      = panel_name,
    n_samples  = nrow(merged_meta),
    n_matched  = matched,
    match_pct  = match_pct,
    n_PD       = sum(merged_meta$phenotype_clean == "PD", na.rm = TRUE),
    n_Control  = sum(merged_meta$phenotype_clean == "Control", na.rm = TRUE),
    n_unmatched = sum(is.na(merged_meta$GP2ID)),
    n_from_manifest = if ("metadata_source" %in% names(merged_meta))
                        sum(merged_meta$metadata_source == "manifest") else 0L,
    n_ancestry_known = sum(!is.na(merged_meta$ancestry)),
    # how many samples are actually usable (with complete covariates).
    n_da_ready = sum(
      merged_meta$phenotype_clean %in% c("PD", "Control") &
      stats::complete.cases(merged_meta[, DA_COVARIATES, drop = FALSE]),
      na.rm = TRUE)
  ))
}


# 3. Merge metadata for each panel 
SE_IN <- list(Neuro220 = SE_FILTERED)
summaries <- list()
for (p in PANELS) summaries[[p]] <- merge_metadata(SE_IN[[p]], p)


# 4. Save merge summary
summary_df <- bind_rows(summaries)
write.csv(summary_df,
          file.path(RESULTS_02_DIR, paste0(COHORT, "_merge_summary.csv")),
          row.names = FALSE)
print(summary_df)
