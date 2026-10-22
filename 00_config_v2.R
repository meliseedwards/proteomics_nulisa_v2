
# Config for NULISA proteomics pipeline (v2)
# Date: October 2026
# Description: Paths & parameters for NULISA GP2 Pilot proteomics workflow




# --- VM Paths -----------------------------------------------------------------
BUCKET_DIR   <- "~/proteomics_bucket/proteomics"
DATA_DIR     <- file.path(BUCKET_DIR, "metadata")
MANIFEST_DIR <- file.path(DATA_DIR, "GP2_nulisa_manifests")
PROJECT_DIR <- "/mnt/output/proteomics/results/nulisa_v2"



# --- Cohort to analyze --------------------------------------------------------
# change cohort and plate correction, source 00, then run scripts. 
COHORT <- "p_ppmi_v10"  
PLATE_CORRECTION <- "covariate" # none, covariate, sva

# "default" = the NA/zero thresholds below; "none" = no NA or zero filtering.
FILTER_MODE <- "none"

# label whether model was adj for plate or not 
PLATE_SUFFIX <- switch(PLATE_CORRECTION, none = "_plate_none", covariate = "_plate_adj", sva = "_plate_sva")

# skip only for ppmi files 
NPQ_SKIP <- 0

# --- Meta-analysis cohorts ---------------------------------------------------
META_COHORTS <- c("ELPD", "UMKLM_KUL", "NTUH", "TRAPCAF", "CANDAS-SMPD")
META_INPUT   <- "plate_adj" # plate_adj = some adjusted cohorts; plate_none = no plate adj any cohort
LOGFC_THRESHOLD <- 0.20   # 0.14 = 10%, 0.20 = 15%

if (COHORT == "p121") {         # NTUH 

  PANELS             <- c("Neuro220")         
  MATRIX_TYPE        <- "PLASMA"
  SAMPLE_DET_MIN     <- 0.90            
  ANCESTRY_KEEP      <- NULL
  VISIT_KEEP         <- NULL            # no visit column = null
  ID_METHOD          <- "strip_suffix"  # SampleName minus "_s2" = GP2ID
  ID_COLUMN          <- NULL
  GP2_JOIN_KEY       <- "GP2ID"
  STUDY_FILTER       <- c("NTUH")
  EXCLUDE_FROM_DA <- c("APOE4")   
  USE_CURATED        <- FALSE
  SAMPLE_ID_PATTERN  <- "^NTUH_[0-9]+_s[0-9]+$"
  MANIFEST_FILES     <- NULL
  MANIFEST_SAMPLE_COL <- NULL
  MANIFEST_GP2ID_COL  <- NULL
  NPQ_DIR  <- file.path(BUCKET_DIR, "P121 GP2")
  NPQ_FILE <- "P121_BSHRI_NULISAseq_Neuro220_NPQ_06092026.xlsx"

} else if (COHORT == "p118") {          # UMKLM AND KUL COMBINED
  PANELS             <- c("Neuro220")
  MATRIX_TYPE        <- "PLASMA"
  SAMPLE_DET_MIN     <- 0.90
  ANCESTRY_KEEP      <- NULL
  VISIT_KEEP         <- NULL
  ID_METHOD          <- "manifest"
  ID_COLUMN          <- NULL
  GP2_JOIN_KEY       <- "GP2ID"
  STUDY_FILTER       <- c("UMKLM", "KUL")
  EXCLUDE_FROM_DA    <- c("APOE4")
  USE_CURATED        <- FALSE
  SAMPLE_ID_PATTERN  <- NULL
  NPQ_DIR            <- file.path(BUCKET_DIR, "P118 GP2")
  NPQ_FILE           <- "P118_BSHRI_NULISAseq_Neuro220_NPQ_06082026.xlsx"
  MANIFEST_FILES     <- c("UMKLM_selfQCV2_2026-01-23_m11.csv",
                          "KUL_selfQCV2_2026-01-23_m4.csv")
  MANIFEST_SAMPLE_COL <- "sample_id"
  MANIFEST_GP2ID_COL  <- "GP2ID"

} else if (COHORT == "p149") {          # CANDAS-SMPD    
  PANELS             <- c("Neuro220")
  MATRIX_TYPE        <- "PLASMA"
  SAMPLE_DET_MIN     <- 0.90
  ANCESTRY_KEEP      <- NULL
  VISIT_KEEP         <- NULL
  ID_METHOD          <- "manifest"
  ID_COLUMN          <- NULL
  GP2_JOIN_KEY       <- "GP2ID"
  STUDY_FILTER       <- c("CANDAS-SMPD")
  EXCLUDE_FROM_DA    <- c("APOE4")
  USE_CURATED        <- FALSE
  SAMPLE_ID_PATTERN  <- NULL
  NPQ_DIR            <- file.path(BUCKET_DIR, "P149 GP2")
  NPQ_FILE           <- "P149_BSHRI_NULISAseq_Neuro220Panel_NPQ_08112026.xlsx"
  MANIFEST_FILES     <- "CANDAS-SMPD_selfQCV2_2026-04-09_m2.csv"
  MANIFEST_SAMPLE_COL <- "sample_id"
  MANIFEST_GP2ID_COL  <- "GP2ID"

} else if (COHORT == "p136") {           # TRAPCAF
  PANELS             <- c("Neuro220")
  MATRIX_TYPE        <- "PLASMA"
  SAMPLE_DET_MIN     <- 0.90
  ANCESTRY_KEEP      <- NULL
  VISIT_KEEP         <- NULL
  ID_METHOD          <- "manifest"
  ID_COLUMN          <- NULL
  GP2_JOIN_KEY       <- "GP2ID"
  STUDY_FILTER       <- c("TRAPCAF")
  EXCLUDE_FROM_DA    <- c("APOE4")
  USE_CURATED        <- FALSE
  SAMPLE_ID_PATTERN  <- NULL
  NPQ_DIR            <- file.path(BUCKET_DIR, "P136 GP2")
  NPQ_FILE           <- "P136_NEW_BSHRI_NULISAseq_Neuro220Panel_NPQ_0825262026.xlsx"
  MANIFEST_FILES     <- c("TRAPCAF_selfQCV2_2026-02-26_m2.csv",
                          "TRAPCAF_selfQCV2_2026-02-26_m3.csv")
  MANIFEST_SAMPLE_COL <- "sample_id"
  MANIFEST_GP2ID_COL  <- "GP2ID"

} else if (COHORT == "p148") {           # PROSPECT
  PANELS             <- c("Neuro220")
  MATRIX_TYPE        <- "PLASMA"
  SAMPLE_DET_MIN     <- 0.90
  ANCESTRY_KEEP      <- NULL
  VISIT_KEEP         <- NULL
  ID_METHOD          <- "manifest"
  ID_COLUMN          <- NULL
  GP2_JOIN_KEY       <- "GP2ID"
  STUDY_FILTER       <- c("PROSPECT")
  EXCLUDE_FROM_DA    <- c("APOE4")
  USE_CURATED        <- FALSE
  SAMPLE_ID_PATTERN  <- NULL
  NPQ_DIR            <- file.path(BUCKET_DIR, "P148 GP2")
  NPQ_FILE           <- "P148_BSHRI_NULISAseq_Neuro220Panel_NPQ_08112026.xlsx"
  MANIFEST_FILES     <- "PROSPECT_selfQCV2_2026-03-27_m19.csv"
  MANIFEST_SAMPLE_COL <- "sample_id"
  MANIFEST_GP2ID_COL  <- "GP2ID"

} else if (COHORT == "p143") {           # Nigeria-PD
  PANELS             <- c("Neuro220")
  MATRIX_TYPE        <- "PLASMA"
  SAMPLE_DET_MIN     <- 0.90
  ANCESTRY_KEEP      <- NULL
  VISIT_KEEP         <- NULL
  ID_METHOD          <- "manifest"
  ID_COLUMN          <- NULL
  GP2_JOIN_KEY       <- "GP2ID"
  STUDY_FILTER       <- c("Nigeria-PD")
  EXCLUDE_FROM_DA    <- c("APOE4")
  USE_CURATED        <- FALSE
  SAMPLE_ID_PATTERN  <- NULL
  NPQ_DIR            <- file.path(BUCKET_DIR, "P143 GP2")
  NPQ_FILE           <- "P143_BSHRI_NULISAseq_Neuro220Panel_NPQ_08112026.xlsx"
  MANIFEST_FILES     <- c("Nigeria-PD_selfQCV2_2026-03-30_m5.csv",
                          "Nigeria-PD_selfQCV2_2026-01-07_m3.csv")
  MANIFEST_SAMPLE_COL <- "sample_id"
  MANIFEST_GP2ID_COL  <- "GP2ID"

} else if (COHORT == "p111") {           # ELPD
  PANELS             <- c("Neuro220")
  MATRIX_TYPE        <- "PLASMA"
  SAMPLE_DET_MIN     <- 0.90
  ANCESTRY_KEEP      <- NULL
  VISIT_KEEP         <- NULL
  ID_METHOD          <- "manifest"
  ID_COLUMN          <- NULL
  GP2_JOIN_KEY       <- "GP2ID"
  STUDY_FILTER       <- c("ELPD")
  EXCLUDE_FROM_DA    <- c("APOE4")
  USE_CURATED        <- FALSE
  SAMPLE_ID_PATTERN  <- NULL
  NPQ_DIR            <- file.path(BUCKET_DIR, "P111 GP2")
  NPQ_FILE           <- "P111_BSHRI_NULISAseq_Neuro220_NPQ_06082026.xlsx"
  MANIFEST_FILES     <- "ELPD_selfQCV2_2026-01-07_m3.csv"
  MANIFEST_SAMPLE_COL <- "GP2sampleID"
  MANIFEST_GP2ID_COL  <- "GP2ID"

} else if (COHORT == "p144") {           # LARGEPD (all controls)  
  PANELS             <- c("Neuro220")
  MATRIX_TYPE        <- "PLASMA"
  SAMPLE_DET_MIN     <- 0.90
  ANCESTRY_KEEP      <- NULL
  VISIT_KEEP         <- NULL
  ID_METHOD          <- "manifest"
  ID_COLUMN          <- NULL
  GP2_JOIN_KEY       <- "GP2ID"
  STUDY_FILTER       <- c("LARGEPD")
  EXCLUDE_FROM_DA    <- c("APOE4")
  USE_CURATED        <- FALSE
  SAMPLE_ID_PATTERN  <- NULL
  NPQ_DIR            <- file.path(BUCKET_DIR, "P144 GP2")
  NPQ_FILE           <- "P144_BSHRI_NULISAseq_Neuro220Panel_NPQ_0814262026.xlsx"
  MANIFEST_FILES     <- "LARGEPD_selfQCV2_2026-04-09_m3.csv"
  MANIFEST_SAMPLE_COL <- "sample_id"
  MANIFEST_GP2ID_COL  <- "GP2ID"

} else if (COHORT == "p146") {           # PDGNRTN (all cases) 
  PANELS             <- c("Neuro220")
  MATRIX_TYPE        <- "PLASMA"
  SAMPLE_DET_MIN     <- 0.90
  ANCESTRY_KEEP      <- NULL
  VISIT_KEEP         <- NULL
  ID_METHOD          <- "manifest"
  ID_COLUMN          <- NULL
  GP2_JOIN_KEY       <- "GP2ID"
  STUDY_FILTER       <- c("PDGNRTN")
  EXCLUDE_FROM_DA    <- c("APOE4")
  USE_CURATED        <- FALSE
  SAMPLE_ID_PATTERN  <- NULL
  NPQ_DIR            <- file.path(BUCKET_DIR, "P146 GP2")
  NPQ_FILE           <- "P146_BSHRI_NULISAseq_Neuro220Panel_NPQ_08122026.xlsx"
  MANIFEST_FILES     <- "PDGNRTN_selfQCV2_2026-04-07_m23.csv"
  MANIFEST_SAMPLE_COL <- "sample_id"
  MANIFEST_GP2ID_COL  <- "GP2ID"

} else if (COHORT == "p144_p146") {   # LARGEPD controls + PDGNRTN cases for combined DA analysis
  PANELS             <- c("Neuro220")
  MATRIX_TYPE        <- "PLASMA"
  SAMPLE_DET_MIN     <- 0.90
  ANCESTRY_KEEP      <- NULL
  VISIT_KEEP         <- NULL
  ID_METHOD          <- "manifest"
  ID_COLUMN          <- NULL
  GP2_JOIN_KEY       <- "GP2ID"
  STUDY_FILTER       <- c("LARGEPD", "PDGNRTN")
  EXCLUDE_FROM_DA    <- c("APOE4")
  USE_CURATED        <- FALSE
  SAMPLE_ID_PATTERN  <- NULL
  NPQ_DIR            <- BUCKET_DIR
  NPQ_FILE           <- c("P144 GP2/P144_BSHRI_NULISAseq_Neuro220Panel_NPQ_0814262026.xlsx",
                      "P146 GP2/P146_BSHRI_NULISAseq_Neuro220Panel_NPQ_08122026.xlsx")
  MANIFEST_FILES     <- c("LARGEPD_selfQCV2_2026-04-09_m3.csv", "PDGNRTN_selfQCV2_2026-04-07_m23.csv")
  MANIFEST_SAMPLE_COL <- c("sample_id", "alternative_id1") 
  MANIFEST_GP2ID_COL  <- "GP2ID"
 
} else if (COHORT == "p_ppmi") {          # ppmi plasma Neuro220 (Project 312), baseline only
  PANELS             <- c("Neuro220")
  MATRIX_TYPE        <- "PLASMA"
  SAMPLE_DET_MIN     <- 0.90
  ANCESTRY_KEEP      <- NULL
  VISIT_KEEP         <- "BL"
  ID_METHOD          <- "column"
  ID_COLUMN          <- "PatNo"
  GP2_JOIN_KEY       <- "clinical_id"     # patno is clinical_id in r12
  STUDY_FILTER       <- c("PPMI-G", "PPMI-N")
  EXCLUDE_FROM_DA    <- c("APOE4", "mCherry") # mCherry is IC (control)
  USE_CURATED        <- FALSE
  SAMPLE_ID_PATTERN  <- NULL
  MANIFEST_FILES     <- NULL
  MANIFEST_SAMPLE_COL <- NULL
  MANIFEST_GP2ID_COL  <- NULL
  NPQ_DIR  <- file.path(BUCKET_DIR, "PPMI", "extracted")
  NPQ_FILE <- "PPMI_Project_312_NULISAseq_Neuro_Panel_1-NPQ_Counts_Plasma_Cohort.xlsx"
  NPQ_SKIP <- 1                            # header note above the column row
} else if (COHORT == "p_ppmi_v10") {       # v10 = follow up time point
  PANELS             <- c("Neuro220")
  MATRIX_TYPE        <- "PLASMA"
  SAMPLE_DET_MIN     <- 0.90
  ANCESTRY_KEEP      <- NULL
  VISIT_KEEP <- "V10"
  ID_METHOD          <- "column"
  ID_COLUMN          <- "PatNo"
  GP2_JOIN_KEY       <- "clinical_id"     
  STUDY_FILTER       <- c("PPMI-G", "PPMI-N")
  EXCLUDE_FROM_DA    <- c("APOE4", "mCherry") 
  USE_CURATED        <- FALSE
  SAMPLE_ID_PATTERN  <- NULL
  MANIFEST_FILES     <- NULL
  MANIFEST_SAMPLE_COL <- NULL
  MANIFEST_GP2ID_COL  <- NULL
  NPQ_DIR  <- file.path(BUCKET_DIR, "PPMI", "extracted")
  NPQ_FILE <- "PPMI_Project_312_NULISAseq_Neuro_Panel_1-NPQ_Counts_Plasma_Cohort.xlsx"
  NPQ_SKIP <- 1                           
} else {
  stop("Unknown COHORT: ", COHORT)
}

# Results directories 
COHORT_DIR <- switch(COHORT,
  p111  = "results_p111_v2",
  p118  = "results_p118_v2",
  p121  = "results_p121_v2",
  p136  = "results_p136_v2",
  p143  = "results_p143_v2",
  p148  = "results_p148_v2",
  p149  = "results_p149_v2",
  p144  = "results_p144_v2",
  p146  = "results_p146_v2",
  p144_p146 = "results_p144_p146_v2",
  p_ppmi = "results_p_ppmi_v2",
  p_ppmi_v10 = "results_p_ppmi_v10_v2",
  stop("No COHORT_DIR defined for ", COHORT))

RESULTS_BASE <- file.path(PROJECT_DIR, paste0(COHORT_DIR, "_filter_", FILTER_MODE))
RESULTS_DIR  <- paste0(RESULTS_BASE, PLATE_SUFFIX) # used for carrier analyses only

RESULTS_01_DIR <- file.path(RESULTS_BASE, "01_load_filter")
RESULTS_02_DIR <- file.path(RESULTS_BASE, "02_merge_gp2")
RESULTS_03_DIR <- file.path(RESULTS_BASE, paste0("03_differential_abundance", PLATE_SUFFIX))
RESULTS_04_DIR <- file.path(RESULTS_BASE, paste0("04_qc_report", PLATE_SUFFIX))

# SC drift script doesn't depend on the DA model, so it lives in RESULTS_BASE
SC_DRIFT_DIR  <- file.path(RESULTS_BASE, "02b_sc_drift")
SC_DRIFT_FILE <- file.path(SC_DRIFT_DIR, paste0("sc_drift_", COHORT, ".csv"))


# --- NULISA data files ----------------------------------------------------------
# release 12
MASTER_KEY_FILE <- "INTERNAL_USE_ONLY_master_key_release12_final_vwb.csv"

# Carrier information in variant report
VARIANT_REPORT_FILE <- "variant_report_r12_final_aug3.csv"


# --- Meta6 GWAS gene list ---------------------------------------------------
META6_FILE <- file.path(DATA_DIR, "meta6_genes.csv")



# --- Pipeline output files (downstream scripts reference these) -------------

# The merged SE for whichever cohort is active. Scripts 04, 06, 08 read this.
SE_FILTERED <- file.path(RESULTS_01_DIR, paste0("se_", COHORT, "_Neuro220_filtered.rds"))
SE_MERGED   <- file.path(RESULTS_02_DIR, paste0("se_", COHORT, "_Neuro220_with_metadata.rds"))


# --- Shared parameters ------------------------------------------------------
FDR_THRESHOLD <- 0.05
OUTLIER_SD_THRESHOLD <- 4


# --- v2 target filters ------------------------------------------------------
ZERO_PROP_MAX <- if (FILTER_MODE == "none") 1 else 0.55
NA_PROP_MAX   <- if (FILTER_MODE == "none") 1 else 0.70

# Does a zero count toward the NA total?
#   FALSE = zeros and NAs are separate states (our reading of the NPQ formula)
ZEROS_COUNT_AS_NA <- FALSE

# Target-level detectability (>=50%), matrix-independent, per Alamar.
# v2 CHANGE: this is now a flag carried into the DA results, not a filter.
DET_THRESH <- 0.50
APPLY_DET_FILTER <- FALSE

# Samples below this detectability (%) ARE excluded from DA and QC. 
# Note Alamar's floor is stricter and matrix-dependent (90% plasma, 70% CSF)
SAMPLE_DETECT_EXCLUDE <- 20

# Minimum non-missing samples before a target is modelled. 
MIN_SAMPLES_PER_TARGET <- 20 # none in our data, placeholder

# Track target-detectability and flag those under 50%
REPORT_MIN_DETECT <- 0.50

# --- v2 sample handling -----------------------------------------------------
# v2 CHANGE: SampleQC == "warning" samples are KEPT and flagged
DROP_SAMPLE_QC_WARN <- FALSE
REMOVE_OUTLIERS <- FALSE

# Some donors contributed more than one sample. Until we know what the visit
# labels mean (0W / 52W / unlabelled) and how to handle donors enrolled under 
# multiple sub-study IDs, these donors are excluded rather than picked between.
DROP_REPEAT_DONORS <- TRUE

# Samples where the proteomic APOE4 readout disagrees with the donor's APOE
# genotype, indicating a probable sample swap.
IDENTITY_EXCLUDE_FILE <- "~/proteomics_nulisa/scripts/identity_flagged_samples.csv"
IDENTITY_EXCLUDE <- if (file.exists(IDENTITY_EXCLUDE_FILE)) {
  x <- read.csv(IDENTITY_EXCLUDE_FILE, stringsAsFactors = FALSE)
  x$SampleName[x$APOE4_SAMPLE_CHECK == "FAIL"]
} else character(0)


# --- v2 DA model ------------------------------------------------------------
# Covariates: age and sex only, per Guglielmo. PDTRTMNT is deliberately absent.
DA_COVARIATES <- c("age", "sex_clean")
EXCLUDE_DIAGNOSES <- c("Prodromal", "Prodromal motor PD", "Prodromal non-motor PD",
                       "Unaffected with hyposmia")



# --- Special target categories (per NULISAseqR user guide) -------------------
EXEMPT_RARE_CASE <- TRUE

HIGH_ABUNDANCE_TARGETS <- list(
  Inflammation = c("CRP", "KNG1"),
  CNS_Disease  = c("APOE", "CRP"),
  Neuro220     = c("APOA1", "APOA2", "APOE", "APOH", "C1q",
                   "CRP", "PPBP", "PROS1", "SERPINA3")
)

RARE_CASE_TARGETS <- list(
  Inflammation = character(0),
  CNS_Disease  = character(0),
  Neuro220     = c("APOE4", "mHTT-exon1", "pLRRK2-S1292",
                   "pPRKN-S65", "pQ-ATXN3", "pQ-HTT", "pRAB10-T73",
                   "pRAB12-S106", "pRAB29-T71")
)

# --- Vendor QC flags (Alamar Neuro 220 data sheet, Table 3) ------------------
# All of these are flags will be carried into DA results. Not used as filters.

# Cross-reactivity > 1%. Named vector = the reported XR percentage.
CROSS_REACTIVE_TARGETS <- list(
  Neuro220 = c("CCL4" = 1438.1, "CCL3" = 294.3, "Oligo-SNCA" = 99.3,
               "PGK1" = 22.7, "pPRKN-S65" = 10.1, "VSNL1" = 9.9,
               "pQ-ATXN3" = 5.1, "pQ-HTT" = 1.8),
  CNS_Disease  = character(0),
  Inflammation = character(0)
)

# Subset whose cross-reactive partner is ALSO on the panel.
XR_PARTNER_ON_PANEL <- list(
  Neuro220 = c("Oligo-SNCA", "pPRKN-S65", "pQ-ATXN3", "pQ-HTT")
)

# Vendor plasma detectability < 50% in NORMAL controls, excluding rare-case.
VENDOR_LOW_DETECT <- list(
  Neuro220 = c("YWHAZ", "PTN", "SNCB", "SNAP25", "Aβ38", "CGRP",
               "FCGR2B", "UCHL1", "GSDMD", "PARP1", "UBB")
)

# Vendor plasma CV > 20% (intra or inter) per Alamar Data Sheet (Neuro220).
VENDOR_HIGH_CV_PLASMA <- list(
  Neuro220 = c("pPRKN-S65", "pRAB29-T71", "pQ-HTT", "pRAB12-S106",
               "BASP1", "CGRP", "SNCB", "CST5", "pLRRK2-S1292",
               "pRAB10-T73", "pQ-ATXN3", "PARK7")
)

# Hemoglobin interference (Alamar Inflammation 250 validation; same targets assumed
# to behave the same on Neuro220 - not vendor-confirmed). Value = lowest tested Hb
# concentration (g/L) at which the assay deviated >30%. Flag, not filter.
VENDOR_HEMOLYSIS_SENSITIVE <- list(
  Neuro220 = c("CCL5" = 2, "FLT1" = 2, "IL18" = 2, "S100A12" = 2, "IL16" = 10)
)


# guards
stopifnot(
  exists("PANELS"), exists("MATRIX_TYPE"),
  exists("SAMPLE_ID_PATTERN"), exists("VISIT_KEEP"), exists("ID_METHOD"),
  exists("NPQ_DIR"), exists("NPQ_FILE"),
  exists("COHORT_DIR"), !is.null(COHORT_DIR),
  ID_METHOD %in% c("column", "strip_suffix", "manifest"),
  (ID_METHOD == "manifest") == !is.null(MANIFEST_FILES),
  is.null(SAMPLE_ID_PATTERN) || nzchar(SAMPLE_ID_PATTERN),
  length(SAMPLE_DET_MIN) == 1, SAMPLE_DET_MIN > 0, SAMPLE_DET_MIN <= 1,
  FILTER_MODE %in% c("default", "none"),
  PLATE_CORRECTION %in% c("none", "covariate", "sva")
)

cat("CONFIG | cohort:", COHORT, "| matrix:", MATRIX_TYPE,
    "| panels:", paste(PANELS, collapse = ", "),
    "| plate adj:", PLATE_CORRECTION,
    "| id pattern:", ifelse(is.null(SAMPLE_ID_PATTERN), "none", SAMPLE_ID_PATTERN), "\n")