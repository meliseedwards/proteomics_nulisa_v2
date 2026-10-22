# Script 06b: age and sex interaction analyses
# For each cohort and protein, fit the main model plus an interaction term,
# then meta-analyse the interaction coefficients (random effects, REML) as in 06.

library(tidyverse)
library(SummarizedExperiment)
library(metafor)

source("~/proteomics_nulisa/scripts/nulisa_pipeline_v2/quick_pass/00_config_v2.R")

OUT_DIR <- file.path(PROJECT_DIR, "05b_interaction")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cohorts   <- c(p111 = "ELPD", p118 = "UMKLM_KUL", p121 = "NTUH", p136 = "TRAPCAF", p149 = "CANDAS-SMPD")
use_plate <- c(p111 = TRUE,   p118 = TRUE,        p121 = FALSE,  p136 = TRUE,      p149 = TRUE)

# 1. Per-cohort interaction models
fit_cohort <- function(co) {
  se <- readRDS(file.path(PROJECT_DIR, paste0("results_", co, "_v2_filter_none"), "02_merge_gp2",
                          paste0("se_", co, "_Neuro220_with_metadata.rds")))
  m    <- as.data.frame(colData(se))
  keep <- m$phenotype_clean %in% c("PD", "Control") & !m$SampleName %in% IDENTITY_EXCLUDE &
          !(!is.na(m$pct_above_lod) & m$pct_above_lod < SAMPLE_DETECT_EXCLUDE)
  m    <- m[keep, ]
  x    <- assay(se, "npq")[setdiff(rownames(se), "APOE4"), keep]
  stopifnot(!anyNA(x))

  d <- data.frame(phenotype = factor(m$phenotype_clean, levels = c("Control", "PD")),
                  sex = factor(case_when(tolower(m$sex_clean) %in% c("female", "f") ~ "Female",
                   tolower(m$sex_clean) %in% c("male", "m")   ~ "Male"),
                  levels = c("Female", "Male")),
                  age10     = (m$age - mean(m$age)) / 10, # centered, 10 yrs
                  plate     = factor(m$PlateID))
  stopifnot(!anyNA(d$sex))
  plate_term <- if (use_plate[[co]]) "+ plate" else ""
  f_sex <- as.formula(paste("y ~ phenotype * sex + age10", plate_term))
  f_age <- as.formula(paste("y ~ phenotype * age10 + sex", plate_term))

  bind_rows(lapply(rownames(x), function(tg) {
    d$y <- x[tg, ]
    cs <- summary(lm(f_sex, data = d))$coefficients["phenotypePD:sexMale", ]
    ca <- summary(lm(f_age, data = d))$coefficients["phenotypePD:age10", ]
    tibble(Target = tg, sex_est = cs[1], sex_se = cs[2], age_est = ca[1], age_se = ca[2])
  })) %>% mutate(cohort = cohorts[[co]], n = nrow(d), .before = 1)
}

per_cohort <- bind_rows(lapply(names(cohorts), fit_cohort))
stopifnot(all(table(per_cohort$Target) == length(cohorts)))
print(per_cohort %>% distinct(cohort, n))
write.csv(per_cohort, file.path(OUT_DIR, "interaction_per_cohort.csv"), row.names = FALSE)

# 2. Meta-analysis of each interaction term
meta_term <- function(est, se) {
  per_cohort %>%
    group_by(Target) %>%
    group_modify(~ {
      fit <- rma(yi = .x[[est]], sei = .x[[se]], method = "REML")
      tibble(estimate = as.numeric(fit$b), ci_lb = fit$ci.lb, ci_ub = fit$ci.ub,
             P = fit$pval, I2 = fit$I2, tau = sqrt(fit$tau2))
    }) %>%
    ungroup() %>%
    mutate(FDR = p.adjust(P, method = "BH")) %>%
    arrange(P)
}

meta_sex <- meta_term("sex_est", "sex_se")
meta_age <- meta_term("age_est", "age_se")
write.csv(meta_sex, file.path(OUT_DIR, "meta_interaction_sex.csv"), row.names = FALSE)
write.csv(meta_age, file.path(OUT_DIR, "meta_interaction_age.csv"), row.names = FALSE)

# 3. Report
for (nm in c("sex", "age")) {
  m <- if (nm == "sex") meta_sex else meta_age
  cat("\n", nm, "x phenotype | FDR < 0.05:", sum(m$FDR < FDR_THRESHOLD),
      "| nominal p < 0.05:", sum(m$P < 0.05), "of", nrow(m), "\n")
  print(m %>% mutate(across(c(estimate, ci_lb, ci_ub, tau), ~ round(.x, 3)), I2 = round(I2)) %>% head(10))
}