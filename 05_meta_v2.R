# =============================================================================
# Script 05 (v2): Random-effects meta-analysis of per-cohort DA results
# Date: September 2026
# Description: For each target, combines the cohort-level per-target
#              linear-model estimates (logFC = PD - Control, ordinary SE)
#              with random-effects meta-analysis (metafor::rma, REML tau^2).
# =============================================================================



# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)
library(metafor)

source("~/proteomics_nulisa/scripts/nulisa_pipeline_v2/quick_pass/00_config_v2.R")

stopifnot(META_INPUT %in% c("plate_adj", "plate_none", "sva"))
meta_cohorts <- META_COHORTS

META_DIR <- file.path(PROJECT_DIR, paste0("05_meta_", META_INPUT))
dir.create(META_DIR, showWarnings = FALSE, recursive = TRUE)

cat("META | input:", META_INPUT, "| cohorts:", paste(meta_cohorts, collapse = ", "), "\n")



# --- 1. Load the combined DA table -------------------------------------------

da_all <- read.csv(file.path(PROJECT_DIR, paste0("DA_all_cohorts_combined_", META_INPUT, ".csv")),
                   stringsAsFactors = FALSE)

# Alert if missing cohorts
missing <- setdiff(meta_cohorts, unique(da_all$cohort_label))
if (length(missing) > 0) stop("Cohort(s) not in combined file: ", paste(missing, collapse = ", "))

da <- da_all %>%
  filter(cohort_label %in% meta_cohorts, !is.na(logFC), !is.na(SE), SE > 0)

# One row per target x cohort
dup_check <- da %>% dplyr::count(Target, cohort_label) %>% filter(n > 1)
if (nrow(dup_check) > 0) { print(dup_check); stop("Duplicate Target x cohort rows in input.") }

cat("Rows in:", nrow(da), "| targets:", n_distinct(da$Target), "\n")
print(da %>% dplyr::count(cohort_label, plate_model, name = "n_targets"))
cat("SE range:", signif(range(da$SE), 3), "| logFC range:", signif(range(da$logFC), 3), "\n")



# --- 2. One target -----------------------------------------------------------
# Input = one row per cohort for a single target
meta_one <- function(d) {
  k <- nrow(d)
  if (k < 2) return(tibble(k = k))

  fit_kh <- rma(yi = d$logFC, sei = d$SE, method = "REML", test = "knha")
  fit_z  <- rma(yi = d$logFC, sei = d$SE, method = "REML")
  fit_fe <- rma(yi = d$logFC, sei = d$SE, method = "FE")

  pooled <- as.numeric(fit_kh$b)
  pred   <- predict(fit_kh)
  loo    <- leave1out(fit_kh)
  shift  <- abs(loo$estimate - pooled)

  tibble(k         = k,
         logFC     = pooled,
         ci_lb     = fit_kh$ci.lb,            # KH (conservative) CI
         ci_ub     = fit_kh$ci.ub,
         P_KH      = fit_kh$pval,
         P_z       = fit_z$pval,
         P_FE      = fit_fe$pval,
         I2        = fit_kh$I2,
         tau2      = fit_kh$tau2,
         pi_lb     = pred$pi.lb,              # prediction interval
         pi_ub     = pred$pi.ub,
         loo_cohort    = d$cohort_label[which.max(shift)],
         loo_same_sign = all(sign(loo$estimate) == sign(pooled)),
         n_sig     = sum(d$adj.P.Val < FDR_THRESHOLD, na.rm = TRUE),
         n_up      = sum(d$logFC > 0),
         n_down    = sum(d$logFC < 0))
}



# --- 3. Run every target -----------------------------------------------------

meta <- da %>%
  group_by(Target) %>%
  group_modify(~ meta_one(.x)) %>%
  ungroup() %>%
  mutate(FDR_KH = p.adjust(P_KH, method = "BH"),
         FDR_z  = p.adjust(P_z,  method = "BH"),
         FDR_FE = p.adjust(P_FE, method = "BH"),
         direction = case_when(FDR_z < FDR_THRESHOLD & logFC > 0 ~ "Up in PD",
                               FDR_z < FDR_THRESHOLD & logFC < 0 ~ "Down in PD",
                               TRUE ~ "Not significant")) %>%
  arrange(P_z)

# Per-target QC context from the cohort files
qc_ctx <- da %>%
  group_by(Target) %>%
  summarise(min_detectability  = if (all(is.na(detectability))) NA_real_ else min(detectability, na.rm = TRUE),
            any_low_detect     = any(low_detect_flag, na.rm = TRUE),
            any_vendor_high_cv = any(vendor_high_cv, na.rm = TRUE),
            any_cross_reactive = any(cross_reactive, na.rm = TRUE),
            any_zero_imbalance = any(zero_imbalance_flag, na.rm = TRUE),
            any_hemolysis      = any(hemolysis_sensitive, na.rm = TRUE),
            in_meta6           = any(in_meta6),
            reportable_all     = all(reportable),
            .groups = "drop")

meta <- meta %>% left_join(qc_ctx, by = "Target") %>% mutate(meta_input = META_INPUT)
stopifnot(!any(duplicated(meta$Target)))

write.csv(meta, file.path(META_DIR, paste0("meta_results_", META_INPUT, ".csv")), row.names = FALSE)



# --- 4. Report ---------------------------------------------------------------

cat("\nTargets meta-analysed:", nrow(meta), "| cohorts per target:", paste(names(table(meta$k)), collapse = ","), "\n")
cat("FDR <", FDR_THRESHOLD, "| Knapp-Hartung:", sum(meta$FDR_KH < FDR_THRESHOLD, na.rm = TRUE),
    "| conventional RE:", sum(meta$FDR_z < FDR_THRESHOLD, na.rm = TRUE),
    "| fixed effect:", sum(meta$FDR_FE < FDR_THRESHOLD, na.rm = TRUE), "\n")
cat("Median I2:", round(median(meta$I2, na.rm = TRUE), 1),
    "| targets with I2 > 50%:", sum(meta$I2 > 50, na.rm = TRUE), "\n")
cat("Conventional-RE hits that keep their sign in every leave-one-out refit:",
    sum(meta$FDR_z < FDR_THRESHOLD & meta$loo_same_sign, na.rm = TRUE), "\n")
cat("Conventional-RE hits reportable in every cohort (detectability >=", REPORT_MIN_DETECT, "):",
    sum(meta$FDR_z < FDR_THRESHOLD & meta$reportable_all, na.rm = TRUE), "\n")

cat("\nTop 20 (by conventional RE p):\n")
print(meta %>%
        select(Target, k, logFC, ci_lb, ci_ub, P_KH, FDR_KH, P_z, FDR_z, P_FE, I2, pi_lb, pi_ub,
               n_sig, n_up, n_down, loo_same_sign, loo_cohort, reportable_all) %>%
        mutate(across(c(logFC, ci_lb, ci_ub, pi_lb, pi_ub), ~ round(.x, 2)), I2 = round(I2, 0)) %>%
        head(20), n = 20, width = Inf)



# --- 5. Figures --------------------------------------------------------------

# forest plots: per-cohort estimates + pooled diamond (KH CI)
top <- meta %>% filter(reportable_all) %>% arrange(P_z) %>% head(12) %>% pull(Target)
png(file.path(META_DIR, paste0("forest_top_hits_", META_INPUT, ".png")),
    width = 14, height = 10, units = "in", res = 200)
par(mfrow = c(3, 4), mar = c(4, 4, 3, 1))
for (tg in top) {
  d <- da %>% filter(Target == tg)
  fit <- rma(yi = d$logFC, sei = d$SE, method = "REML", test = "knha")
  forest(fit, slab = d$cohort_label, xlab = "log2FC (PD vs Control)", main = tg, cex = 0.9)
}
dev.off()

# meta volcano
p_vol <- ggplot(meta, aes(logFC, -log10(P_z), colour = direction)) +
  geom_point(alpha = 0.7, size = 2.2) +
  geom_point(data = filter(meta, I2 > 50), shape = 1, size = 4, colour = "grey30") +
  ggrepel::geom_text_repel(data = filter(meta, FDR_z < FDR_THRESHOLD, reportable_all) %>% head(20),
                           aes(label = Target), colour = "black", size = 3, max.overlaps = 30) +
  scale_colour_manual(values = c("Up in PD" = "tomato", "Down in PD" = "steelblue",
                                 "Not significant" = "grey80")) +
  labs(title = paste0("Meta-analysis (", META_INPUT, "): PD vs Control, ", length(meta_cohorts), " cohorts"),
       subtitle = "random effects, conventional inference; grey ring = I2 > 50%",
       x = "pooled log2FC", y = "-log10(p)") +
  theme_bw()
ggsave(file.path(META_DIR, paste0("volcano_meta_", META_INPUT, ".png")), p_vol,
       width = 9, height = 7, dpi = 150)

# per-cohort log2FC heatmap for the top 30 targets
sig <- meta %>% filter(reportable_all) %>% arrange(P_z) %>% head(30) %>% pull(Target)
hm <- da %>% filter(Target %in% sig) %>%
  mutate(Target = factor(Target, levels = rev(sig)),
         sig_in_cohort = adj.P.Val < FDR_THRESHOLD)
p_hm <- ggplot(hm, aes(cohort_label, Target, fill = logFC)) +
  geom_tile() +
  geom_point(data = filter(hm, sig_in_cohort), shape = 8, size = 1.5) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato", midpoint = 0, name = "log2FC") +
  labs(title = "Top 30 meta-analysis targets: per-cohort log2FC (* = FDR < 0.05 in that cohort)",
       x = NULL, y = NULL) +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(META_DIR, paste0("heatmap_meta_top30_", META_INPUT, ".png")), p_hm,
       width = 8, height = max(4, 0.25 * length(sig) + 2), dpi = 150)

cat("\nWrote results and figures to", META_DIR, "\n")