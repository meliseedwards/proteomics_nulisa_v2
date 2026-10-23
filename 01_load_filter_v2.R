# Script 01 (v2): load, filter, & build SummarizedExperiment (SE) objects
# Date: October 2026
# Description: Loads NPQ data, filters to study samples, builds SE. 



# 0. Setup
# Note: SummarizedExperiment must be loaded before tidyverse
library(SummarizedExperiment)
library(tidyverse)
library(readxl)

source("~/proteomics_nulisa/scripts/nulisa_pipeline_v2/quick_pass/00_config_v2.R")

dir.create(RESULTS_01_DIR, showWarnings = FALSE, recursive = TRUE)


# 1. Process one panel
process_panel <- function(npq_file, panel_name) {
  npq_path <- file.path(NPQ_DIR, npq_file)

  # Special-target categories for this panel 
  ha_targets   <- HIGH_ABUNDANCE_TARGETS[[panel_name]]
  rare_targets <- RARE_CASE_TARGETS[[panel_name]]
  stopifnot(!is.null(ha_targets), !is.null(rare_targets))
  qc_exclude_targets <- c(ha_targets, rare_targets)

  # Load NPQ data
  npq_long <- bind_rows(lapply(npq_path, function(p)
  read_excel(p, sheet = 1, skip = NPQ_SKIP, na = "NA", col_types = "text")))
  cat("Raw rows in the data:", nrow(npq_long), "\n")

  # Repair garbled Aβ target names (PPMI export)
  before <- unique(npq_long$Target)
  npq_long$Target <- fix_target_names(npq_long$Target)
  changed <- setdiff(unique(npq_long$Target), before)
  if (length(changed) > 0)
    cat("fixed", length(changed), "target name(s):", paste(changed, collapse = ", "), "\n")

 # Sites name the same columns differently so rename to one standard set
  if ("SampleMatrixType" %in% names(npq_long))
    npq_long <- dplyr::rename(npq_long, Biofluid = SampleMatrixType)
  if ("targetLOD_NPQ" %in% names(npq_long))
    npq_long <- dplyr::rename(npq_long, LOD = targetLOD_NPQ)
      if ("ClinicalEvent" %in% names(npq_long))
    npq_long <- dplyr::rename(npq_long, CLINICAL_EVENT = ClinicalEvent)


# Donor ID: matching IDs for different cohorts
  # Three routes, set by ID_METHOD in the config:
  #   "column"        the file already names the participant (PPMI: PATNO)
  #   "strip_suffix"  SampleName minus the "_s2" aliquot suffix is the GP2ID
  #   "manifest"      SampleName is a site-local ID; the manifest maps it to
  #                   a GP2ID (P118 holds both UMKLM and KUL)
  if (ID_METHOD == "column") {
    npq_long$DONOR_ID <- as.character(npq_long[[ID_COLUMN]])

  } else if (ID_METHOD == "strip_suffix") {
    npq_long$DONOR_ID <- sub("_s[0-9]+$", "", npq_long$SampleName)

  } else if (ID_METHOD == "manifest") {
    read_one <- function(f) {
      if (grepl("\\.xlsx$", f))
        readxl::read_excel(file.path(MANIFEST_DIR, f), sheet = 2,
                           .name_repair = "unique_quiet")
      else
        read.csv(file.path(MANIFEST_DIR, f), stringsAsFactors = FALSE)
    }
    
    lookup <- bind_rows(lapply(seq_along(MANIFEST_FILES), function(i) {
    m   <- read_one(MANIFEST_FILES[i])
    col <- if (length(MANIFEST_SAMPLE_COL) > 1) MANIFEST_SAMPLE_COL[i] else MANIFEST_SAMPLE_COL
    m %>% transmute(SampleName = trimws(gsub("'", "", as.character(.data[[col]]))),
                      DONOR_ID   = trimws(as.character(.data[[MANIFEST_GP2ID_COL]])))
    })) %>%
      filter(SampleName != "", DONOR_ID != "") %>%
      distinct()

    conflicts <- lookup %>% dplyr::count(SampleName) %>% filter(n > 1)
     if (nrow(conflicts) > 0) {
      print(lookup %>% filter(SampleName %in% conflicts$SampleName) %>% arrange(SampleName))
      stop(nrow(conflicts), " sampleName(s) map to more than one donor ID in the manifests")
     }


    n_before <- n_distinct(npq_long$SampleName)
    npq_long <- npq_long %>%
      mutate(SampleName = trimws(SampleName)) %>%
      left_join(lookup, by = "SampleName")
  
    cat("Manifest matched", n_distinct(npq_long$SampleName[!is.na(npq_long$DONOR_ID)]),
        "of", n_before, "wells\n")
    npq_long <- npq_long %>% filter(!is.na(DONOR_ID))

  } else stop("Unknown ID_METHOD: ", ID_METHOD)

# filter
# toupper() on both sides: PPMI CSF says "CSF"/"passed", NTUH says
# "PLASMA"/"Passed", PPMI plasma says "plasma"/"PASS".
npq_samples <- npq_long %>%
  filter(toupper(SampleType) == "SAMPLE",
     toupper(Biofluid) == toupper(MATRIX_TYPE),
        !is.na(DONOR_ID), DONOR_ID != "NA")

  if (!is.null(VISIT_KEEP)) {
    npq_samples <- npq_samples %>% filter(CLINICAL_EVENT %in% VISIT_KEEP)
  }

  if (!is.null(SAMPLE_ID_PATTERN)) {
    bad_id <- !grepl(SAMPLE_ID_PATTERN, npq_samples$SampleName)
    if (any(bad_id)) {
      cat("Dropping", n_distinct(npq_samples$SampleName[bad_id]),
          "well(s) whose name does not match the cohort ID pattern:",
          paste(unique(npq_samples$SampleName[bad_id]), collapse = ", "), "\n")
      npq_samples <- npq_samples[!bad_id, ]
    }
  }
  cat("After filtering to study samples:",
      length(unique(npq_samples$SampleName)), "samples\n")

  # report the SampleQC split, then only drop if the config says to
  qc_split <- npq_samples %>% distinct(SampleName, SampleQC) %>% dplyr::count(SampleQC)
  cat("SampleQC:", paste(qc_split$SampleQC, qc_split$n, collapse = ", "), "\n")

  if (DROP_SAMPLE_QC_WARN) {
    npq_samples <- npq_samples %>% filter(tolower(SampleQC) %in% c("passed", "pass"))
    cat("  DROP_SAMPLE_QC_WARN = TRUE, so warned samples were removed:",
        length(unique(npq_samples$SampleName)), "samples remain\n")
  }


  # Donors with more than one sample. 
  # do not know what the visit labels mean, so drop the donor for now and record the list.
  repeat_donors <- npq_samples %>%
    distinct(SampleName, DONOR_ID) %>%
    dplyr::count(DONOR_ID) %>%
    dplyr::filter(n > 1) %>%
    pull(DONOR_ID)

  if (length(repeat_donors) > 0) {
    cat("Donors with more than one sample:", length(repeat_donors), "\n")
    write.csv(npq_samples %>%
                filter(DONOR_ID %in% repeat_donors) %>%
                distinct(SampleName, DONOR_ID, PlateID),
              file.path(RESULTS_01_DIR, paste0("repeat_donors_", panel_name, ".csv")),
              row.names = FALSE)
    if (DROP_REPEAT_DONORS) {
      npq_samples <- npq_samples %>% filter(!DONOR_ID %in% repeat_donors)
      cat("  DROP_REPEAT_DONORS = TRUE, so all their samples were removed\n")
    } else {
      stop("DROP_REPEAT_DONORS is FALSE but no selection rule is defined.")
    }
  }

  cat("After repeat handling:", n_distinct(npq_samples$SampleName), "samples |",
      n_distinct(npq_samples$DONOR_ID), "donors\n")
  stopifnot(!any(duplicated(unique(npq_samples[c("SampleName","DONOR_ID")])$DONOR_ID)))
  
  # Coerce NPQ/LOD to numeric
  npq_samples <- npq_samples %>%
    mutate(NPQ = as.numeric(NPQ), LOD = as.numeric(LOD))


  # LOD diagnostic: % of samples below the plate LOD, per target
  # Alamar does not report LOD or detectability for high-abundance targets
  # so their LOD is NA and "below LOD" is undefined
  lod_report <- npq_samples %>%
    mutate(high_abundance = Target %in% ha_targets,
           rare_case      = Target %in% rare_targets) %>%
    group_by(Target, high_abundance, rare_case) %>%
    summarise(pct_below_lod = round(100 * mean(NPQ < LOD, na.rm = TRUE), 1),
              .groups = "drop") %>%
    arrange(desc(pct_below_lod))
  
  write.csv(lod_report,
            file.path(RESULTS_01_DIR, paste0("lod_report_", panel_name, ".csv")),
            row.names = FALSE)
  
  # print >50%-below-LOD targets, EXCLUDING high-abundance
  high_below <- lod_report %>% filter(pct_below_lod > 50, !high_abundance)
  
  cat("Targets below LOD in >50% of samples (excl. high-abundance):",
      nrow(high_below), "\n")
  
  if (nrow(high_below) > 0) print(high_below, n = Inf)

  # Sample-level detectability (FLAG, do not auto-drop) 
  # Per Alamar: sample detectability = % of targets above that sample's LOD.
  # Both high-abundance and rare-case targets are excluded from this metric.
  sample_det <- npq_samples %>%
    filter(!Target %in% qc_exclude_targets) %>%
    group_by(SampleName) %>%
    summarise(pct_above_lod = round(100 * mean(NPQ > LOD, na.rm = TRUE), 1),
              .groups = "drop") %>%
    mutate(detectability_flag = pct_above_lod < SAMPLE_DET_MIN * 100)
  
  n_flag <- sum(sample_det$detectability_flag)

  cat("Samples flagged below", SAMPLE_DET_MIN * 100,
      "% detectability (NOT dropped):", n_flag, "\n")
  
  if (n_flag > 0) {
    print(sample_det %>% filter(detectability_flag) %>% arrange(pct_above_lod))
  }

  write.csv(sample_det,
            file.path(RESULTS_01_DIR, paste0("sample_detectability_", panel_name, ".csv")),
            row.names = FALSE)

  # Ensure one value per SampleName x Target before pivoting wide
  dup_measurements <- npq_samples %>%
    dplyr::count(SampleName, Target) %>%
    filter(n > 1)
  
  if (nrow(dup_measurements) > 0) {
    cat("WARNING: ", nrow(dup_measurements),
        " SampleName x Target pairs have duplicate measurements.\n")
    print(head(dup_measurements))
    stop("Duplicate measurements found - investigate before proceeding.")
  }

  # Pivot from long to wide (proteins as rows, samples as columns)
  npq_wide <- npq_samples %>%
    select(SampleName, Target, NPQ) %>%
    pivot_wider(names_from = SampleName,
                values_from = NPQ)
  cat("Wide matrix:", nrow(npq_wide), "targets x",
      ncol(npq_wide) - 1, "samples\n")

  # Build sample metadata: one row per sample with QC and detectability flags
  sample_meta <- npq_samples %>%
    select(SampleName, PlateID, SampleQC, DONOR_ID, any_of("CLINICAL_EVENT"), 
    Biofluid) %>%
    distinct(SampleName, .keep_all = TRUE) %>%
    left_join(sample_det, by = "SampleName") %>%
    mutate(sample_qc_flag = !tolower(SampleQC) %in% c("passed", "pass")) %>%
    arrange(SampleName)


  # Target inclusion decision 
  # computed from this cohort:
  #   na_prop = proportion of values that are NA (no value reported)
  #   zero_prop = proportion of PRESENT values that are exactly 0 (zero reads)
  #   detectability = proportion of present values above that row's LOD
  target_stats <- npq_samples %>%
    group_by(Target) %>%
    summarise(
      n_values      = dplyr::n(),
      n_present     = sum(!is.na(NPQ)),
      n_zero        = sum(!is.na(NPQ) & NPQ == 0),
      na_prop       = mean(is.na(NPQ)),
      zero_prop     = ifelse(n_present > 0, n_zero / n_present, NA_real_),
      detectability = ifelse(n_present > 0,
                             sum(!is.na(NPQ) & !is.na(LOD) & NPQ > LOD) /
                               sum(!is.na(NPQ) & !is.na(LOD)), NA_real_),
      .groups = "drop"
    )

  if (ZEROS_COUNT_AS_NA) {
    target_stats <- target_stats %>%
      mutate(na_prop = (n_values - n_present + n_zero) / n_values)
    cat("NOTE: ZEROS_COUNT_AS_NA = TRUE, so zeros count toward na_prop.\n")
  }
  
  # Give each target a category and a keep/drop status
  target_status <- target_stats %>%
    mutate(
      category = case_when(
        Target %in% ha_targets   ~ "high_abundance",
        Target %in% rare_targets ~ "rare_case",
        TRUE                     ~ "standard"
      ),
      exempt = EXEMPT_RARE_CASE & category == "rare_case",
      status = case_when(
        n_present == 0                       ~ "dropped_no_data",   # high-abundance
        exempt                               ~ "kept_rare_case",
        na_prop   > NA_PROP_MAX              ~ "dropped_high_na",
        zero_prop > ZERO_PROP_MAX            ~ "dropped_high_zero",
        APPLY_DET_FILTER & !is.na(detectability) &
          detectability < DET_THRESH         ~ "dropped_low_detect",
        TRUE                                 ~ "kept"
      ),
      low_detect_flag = !is.na(detectability) & detectability < DET_THRESH
    )

  targets_pass <- target_status %>%
    filter(status %in% c("kept", "kept_rare_case")) %>%
    pull(Target)


  # report + save the per-target decision 
  cat("Target inclusion decision (", panel_name, "):\n", sep = "")
  print(target_status %>% dplyr::count(category, status))

  dropped <- target_status %>% filter(!Target %in% targets_pass)
  if (nrow(dropped) > 0) {
    cat("Dropped targets:\n")
    print(dropped %>% select(Target, category, status, na_prop, zero_prop, detectability))
  }
  cat("KEPT but under", DET_THRESH * 100, "% detectable (flagged, not dropped):",
      sum(target_status$Target %in% targets_pass & target_status$low_detect_flag), "\n")


  # sanity check: a no-data target should be a known high-abundance target
  unexpected <- target_status %>%
    filter(status == "dropped_no_data", category != "high_abundance")
  if (nrow(unexpected) > 0)
    cat("WARNING: target(s) with no NPQ data that are NOT high-abundance:",
        paste(unexpected$Target, collapse = ", "), "\n")
  
  write.csv(target_status,
            file.path(RESULTS_01_DIR, paste0("target_status_", panel_name, ".csv")),
            row.names = FALSE)
  
  
  # Build the SummarizedExperiment (filtered to included targets)
  se <- build_se(npq_wide, sample_meta, target_status, target_subset = targets_pass)
  cat("SE built:", nrow(se), "targets x", ncol(se), "samples\n\n")
  out_path <- file.path(RESULTS_01_DIR, paste0("se_", COHORT, "_", panel_name, "_filtered.rds"))
  saveRDS(se, out_path)
  return(se)
}

# 2. Build a SummarizedExperiment from wide data
# target QC table is stored in rowData so 03 can use it
build_se <- function(wide_df, sample_meta, target_status, target_subset = NULL) {
  df <- wide_df
  if (!is.null(target_subset)) {
    df <- df %>% filter(Target %in% target_subset)
  }
  expr_matrix <- as.matrix(df %>% select(-Target))
  rownames(expr_matrix) <- df$Target
  sample_meta_aligned <- as.data.frame(sample_meta[
  match(colnames(expr_matrix), sample_meta$SampleName), ])
  rownames(sample_meta_aligned) <- sample_meta_aligned$SampleName
  row_data <- as.data.frame(target_status[
  match(rownames(expr_matrix), target_status$Target), ])
  rownames(row_data) <- row_data$Target
  se <- SummarizedExperiment(
    assays  = list(npq = expr_matrix),
    colData = sample_meta_aligned,
    rowData = row_data
  )
  # check samples and targets line up with their metadata
  stopifnot(
    identical(colnames(se), colData(se)$SampleName),
    identical(rownames(se), rowData(se)$Target),
    !anyNA(colData(se)$DONOR_ID)
  )
  se
}

# 3. Process panels 
NPQ_FILES <- list(Neuro220 = NPQ_FILE)

se_list <- list()
for (p in PANELS) {
  se_list[[p]] <- process_panel(NPQ_FILES[[p]], p)
}

summary_df <- data.frame(
  panel            = PANELS,
  n_samples        = sapply(se_list, ncol),
  n_targets_passed = sapply(se_list, nrow),
  n_sampleqc_flag  = sapply(se_list, function(s) sum(colData(s)$sample_qc_flag, na.rm = TRUE)),
  n_sampledet_flag = sapply(se_list, function(s) sum(colData(s)$detectability_flag, na.rm = TRUE)),
  n_lowdetect_kept = sapply(se_list, function(s) sum(rowData(s)$low_detect_flag, na.rm = TRUE))
  )

# 4. Save summary
write.csv(summary_df, file.path(RESULTS_01_DIR, paste0(COHORT, "_load_summary.csv")),
          row.names = FALSE)
print(summary_df)
