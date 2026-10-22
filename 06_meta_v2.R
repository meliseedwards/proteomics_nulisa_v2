# Script 06: Random-effects meta-analysis of per-cohort DA results
# Date: October 2026
# Description: For each target, we combine the cohort-level per-target
#              linear-model estimates (logFC = PD - Control, ordinary SE)
#              with random-effects meta-analysis (metafor::rma, REML).



# 0. Setup
library(tidyverse)
library(metafor)

source("~/proteomics_nulisa/scripts/nulisa_pipeline_v2/quick_pass/00_config_v2.R")

META_DIR <- file.path(PROJECT_DIR, paste0("05_meta_", META_INPUT))
dir.create(META_DIR, showWarnings = FALSE, recursive = TRUE)

cat("META | input:", META_INPUT, "| cohorts:", paste(META_COHORTS, collapse = ", "), "\n")



# 1. Load the combined DA table 
da_all <- read.csv(file.path(PROJECT_DIR, paste0("DA_all_cohorts_combined_", META_INPUT, ".csv")),
                   stringsAsFactors = FALSE)

# Alert if missing cohorts
missing <- setdiff(META_COHORTS, unique(da_all$cohort_label))
if (length(missing) > 0) stop("Cohort(s) not in combined file: ", paste(missing, collapse = ", "))

da <- da_all %>%
  filter(cohort_label %in% META_COHORTS, !is.na(logFC), !is.na(SE), SE > 0)

# One row per target x cohort
dup_check <- da %>% dplyr::count(Target, cohort_label) %>% filter(n > 1)
if (nrow(dup_check) > 0) { print(dup_check); stop("Duplicate Target x cohort rows in input.") }

cat("Rows in:", nrow(da), "| targets:", n_distinct(da$Target), "\n")
print(da %>% dplyr::count(cohort_label, plate_model, name = "n_targets"))


# 2. Meta-analysis for each protein
meta_one <- function(d) {
  fit    <- rma(yi = d$logFC, sei = d$SE, method = "REML")
  fit_fe <- rma(yi = d$logFC, sei = d$SE, method = "FE")
  pooled <- as.numeric(fit$b)
  loo    <- leave1out(fit)

  tibble(k        = nrow(d),
         logFC    = pooled,
         ci_lb    = fit$ci.lb,
         ci_ub    = fit$ci.ub,
         P_z      = fit$pval,
         I2       = fit$I2,
         tau2     = fit$tau2,
         tau      = sqrt(fit$tau2),
         logFC_FE = as.numeric(fit_fe$b),
         P_FE     = fit_fe$pval,
         loo_same_sign = all(sign(loo$estimate) == sign(pooled)),
         n_up     = sum(d$logFC > 0),
         n_down   = sum(d$logFC < 0))
}


# 3. Run every target 
meta <- da %>%
  group_by(Target) %>%
  group_modify(~ meta_one(.x)) %>%
  ungroup() %>%
  mutate(FDR_z        = p.adjust(P_z,  method = "BH"),
         FDR_FE       = p.adjust(P_FE, method = "BH"),
         pct_change   = round((2^logFC - 1) * 100, 1),
         large_effect = abs(logFC) >= LOGFC_THRESHOLD,
         direction = case_when(FDR_z < FDR_THRESHOLD & logFC > 0 ~ "Up in PD",
                               FDR_z < FDR_THRESHOLD & logFC < 0 ~ "Down in PD",
                               TRUE ~ "Not significant")) %>%
  arrange(P_z)

# QC flags carried over from cohort DA files
qc_ctx <- da %>%
  group_by(Target) %>%
  summarise(any_low_detect     = any(low_detect_flag, na.rm = TRUE),
            any_vendor_high_cv = any(vendor_high_cv, na.rm = TRUE),
            any_zero_imbalance = any(zero_imbalance_flag, na.rm = TRUE),
            any_hemolysis      = any(hemolysis_sensitive, na.rm = TRUE),
            reportable_all     = all(!is.na(reportable) & reportable),
            .groups = "drop")

meta <- meta %>% left_join(qc_ctx, by = "Target") %>% mutate(meta_input = META_INPUT)
stopifnot(!any(duplicated(meta$Target)))

write.csv(meta, file.path(META_DIR, paste0("meta_results_", META_INPUT, ".csv")), row.names = FALSE)


# 4. Report 
cat("\nFDR <", FDR_THRESHOLD, "| random effects:", sum(meta$FDR_z < FDR_THRESHOLD),
    "| fixed effect:", sum(meta$FDR_FE < FDR_THRESHOLD),
    "| median I2:", round(median(meta$I2), 1), "\n")
cat("RE hits that keep their sign when any one cohort is dropped:",
    sum(meta$FDR_z < FDR_THRESHOLD & meta$loo_same_sign), "\n")

print(meta %>%
        select(Target, logFC, pct_change, ci_lb, ci_ub, FDR_z, FDR_FE, I2, tau, n_up, n_down, reportable_all) %>%
        mutate(across(c(logFC, ci_lb, ci_ub, tau), ~ round(.x, 2)), I2 = round(I2)) %>%
        head(20), width = Inf)


# 5. Figures 
# forest plots for top targets
top <- meta %>% filter(FDR_z < FDR_THRESHOLD, reportable_all) %>% head(12) %>% pull(Target)
png(file.path(META_DIR, paste0("forest_top_hits_", META_INPUT, ".png")),
    width = 14, height = 10, units = "in", res = 200)
par(mfrow = c(3, 4), mar = c(4, 4, 3, 1))
for (tg in top) {
  d <- da %>% filter(Target == tg)
  fit <- rma(yi = d$logFC, sei = d$SE, method = "REML")
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
  labs(title = paste0("Meta-analysis (", META_INPUT, "): PD vs Control, ", length(META_COHORTS), " cohorts"),
       subtitle = "random effects; grey ring = I2 > 50%",
       x = "log2 fold change (PD vs Control)", y = "-log10(p)") +
  geom_vline(xintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD), linetype = "dashed", colour = "grey50") +
  theme_bw()
ggsave(file.path(META_DIR, paste0("volcano_meta_", META_INPUT, ".png")), p_vol,
       width = 9, height = 7, dpi = 150)

# per-cohort log2FC heatmap for the top 10 targets
sig <- meta %>% filter(reportable_all) %>% arrange(P_z) %>% head(10) %>% pull(Target)
hm <- da %>% filter(Target %in% sig) %>%
  mutate(Target = factor(Target, levels = rev(sig)),
         sig_in_cohort = adj.P.Val < FDR_THRESHOLD)
p_hm <- ggplot(hm, aes(cohort_label, Target, fill = logFC)) +
  geom_tile() +
  geom_point(data = filter(hm, sig_in_cohort), shape = 8, size = 1.5) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato", midpoint = 0, name = "log2FC") +
  labs(title = "10 lowest-p reportable targets: per-cohort log2FC (* = FDR < 0.05 in that cohort)",
       x = NULL, y = NULL) +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(META_DIR, paste0("heatmap_meta_top10_", META_INPUT, ".png")), p_hm,
       width = 8, height = max(4, 0.25 * length(sig) + 2), dpi = 150)