# Script 02b: Sample control (SC) drift across plates, per target
# Date: August 2026
# Description: SC wells are pooled plasma that aren't used in normalization,
#              so plate-to-plate differences in them are technical, not
#              biological. Reports per-plate medians, drift from the
#              across-plate median, and intra/inter-plate CV per target.


# 0. Setup 

library(SummarizedExperiment)
library(tidyverse)
library(readxl)

source("~/proteomics_nulisa/scripts/nulisa_pipeline_v2/quick_pass/00_config_v2.R")

dir.create(SC_DRIFT_DIR, showWarnings = FALSE, recursive = TRUE)


# 1. Load SC wells
npq <- bind_rows(lapply(file.path(NPQ_DIR, NPQ_FILE), function(f) {
  read_excel(f, sheet = 1, skip = NPQ_SKIP, na = "NA") %>%
    mutate(across(any_of(c("NPQ", "targetLOD_NPQ", "LOD")), as.numeric))
})) %>%
  dplyr::rename(any_of(c(LOD = "targetLOD_NPQ")))

npq$Target <- fix_target_names(npq$Target)

sc <- npq %>%
  dplyr::filter(SampleType == "SC") %>%
  mutate(NPQ    = as.numeric(NPQ),
         LOD    = as.numeric(LOD),
         linear = 2^NPQ - 1, # note: CV is computed on unlogged values
         plate  = paste0("p", match(PlateID, sort(unique(PlateID)))))

cat("SC wells:", dplyr::n_distinct(sc$SampleName),
    "across", dplyr::n_distinct(sc$plate), "plates |",
    dplyr::n_distinct(sc$Target), "targets\n")

# Below-LOD values are excluded from CV per Alamar (per-target, not well) 
# LOD is NA for high-abundance targets, so those are kept
ok <- sc %>% dplyr::filter(is.na(LOD) | NPQ > LOD)


# 2. Per-plate median and intra-plate CV
# For each target on each plate, take the median of the 3 SC replicates
# then measure how far each plate sits from the median across all plates.
per_plate <- ok %>%
  group_by(Target, plate) %>%
  summarise(n_sc_above_lod = dplyr::n(),
            npq_median     = median(NPQ),
            npq_sd         = sd(NPQ),
            intra_cv       = if (dplyr::n() >= 2) 100 * sd(linear) / mean(linear)
                             else NA_real_,
            .groups = "drop") %>%
  group_by(Target) %>%
  mutate(delta_npq = npq_median - median(npq_median)) %>%
  ungroup()

intra_summary <- per_plate %>%
  group_by(Target) %>%
  summarise(intra_cv_median = median(intra_cv, na.rm = TRUE),
            intra_cv_max = if (all(is.na(intra_cv))) NA_real_ else max(intra_cv, na.rm = TRUE),
            worst_plate  = if (all(is.na(intra_cv))) NA_character_ else plate[which.max(intra_cv)],
            .groups = "drop")


# 3. Inter-plate CV across all SC wells 
inter <- ok %>%
  group_by(Target) %>%
  summarise(n_sc_total = dplyr::n(),
            n_plates   = dplyr::n_distinct(plate),
            inter_cv   = if (dplyr::n() >= 2) 100 * sd(linear) / mean(linear)
                         else NA_real_,
            .groups = "drop")

drift_max <- per_plate %>%
  group_by(Target) %>%
  summarise(max_abs_drift = max(abs(delta_npq)), .groups = "drop")


# 4. Assemble 
sc_drift <- per_plate %>%
  pivot_wider(names_from  = plate,
              values_from = c(n_sc_above_lod, npq_median, npq_sd,
                              delta_npq, intra_cv)) %>%
  left_join(intra_summary, by = "Target") %>%
  left_join(drift_max,    by = "Target") %>%
  left_join(inter,        by = "Target") %>%
  arrange(desc(max_abs_drift))


# 5. Report in the terminal 
cat("\nMedian intra-plate SC CV by plate (Alamar target: <10%):\n")
print(per_plate %>%
        group_by(plate) %>%
        summarise(median_cv = round(median(intra_cv, na.rm = TRUE), 2),
                  .groups = "drop"))

cat("\nMedian inter-plate SC CV:",
    round(median(sc_drift$inter_cv, na.rm = TRUE), 2), "% (Alamar target: <15%)\n")

cat("Targets with no inter-plate CV (fewer than 2 SC values above LOD):",
    sum(is.na(sc_drift$inter_cv)), "\n")

cat("Targets not measurable above LOD on all plates:",
    sum(sc_drift$n_plates < dplyr::n_distinct(sc$plate), na.rm = TRUE), "\n")

cat("\nLargest absolute drift from the across-plate median:\n")
print(sc_drift %>%
        select(Target, max_abs_drift, inter_cv, n_sc_total, n_plates) %>%
        head(15))

cat("Targets with at least one plate over Alamar's 10% intra-plate CV:",
    sum(sc_drift$intra_cv_max > 10, na.rm = TRUE), "\n")


#  6. Save 
# Long format is easier to join onto DA results; wide is easier to eyeball.
write.csv(per_plate,
          file.path(SC_DRIFT_DIR, paste0("sc_drift_long_", COHORT, ".csv")),
          row.names = FALSE)

write.csv(sc_drift,
          file.path(SC_DRIFT_DIR, paste0("sc_drift_", COHORT, ".csv")),
          row.names = FALSE)

cat("\nWrote sc_drift_", COHORT, ".csv and sc_drift_long_", COHORT, ".csv\n", sep = "")


#  7. Plate labels 
# p1/p2/p3 are display labels built here from sorted PlateID. Print the
# mapping so the columns can be traced back.
cat("\nPlate label mapping:\n")
print(sc %>% distinct(plate, PlateID) %>% arrange(plate))