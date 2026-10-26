# GP2 NULISA proteomics pipeline

Quality control, differential abundance and meta-analysis for the GP2 NULISA pilot: plasma, Alamar NULISAseq Neuro220 panel (~220 protein targets).

Each cohort is analysed separately (Parkinson's disease cases vs controls). The per-target effect sizes and standard errors are then combined across cohorts by meta-analysis. PPMI (Project 312) runs through the same scripts as an independent replication dataset and is not pooled.

## Approach

Differential abundance runs first; quality control annotates the results afterwards. Every target that can be modelled gets an estimate, and QC flags travel with it. A poorly measured target is flagged, not silently dropped.

Removed before modelling: assay control wells, samples with no resolvable donor ID, donors with more than one sample, samples that are not PD or control or lack covariates, probable sample swaps (proteomic APOE4 readout disagrees with APOE genotype), and samples below 20% detectability (total assay failure).

## Scripts

Run in order. Set `COHORT` and `PLATE_CORRECTION` in `00_config_v2.R`, then run 01–04 for that cohort. Repeat for each cohort, then run 05–06b once.

| Script | Purpose |
|---|---|
| `00_config_v2.R` | Cohort settings, file paths, thresholds, target categories, meta-analysis cohorts |
| `01_load_filter_v2.R` | Load NPQ, resolve donor IDs, compute detectability, build a SummarizedExperiment |
| `02a_merge_v2.R` | Join GP2 release metadata (manifest fallback), harmonise phenotype and sex |
| `02b_SC_drift_v2.R` | Plate drift and intra-/inter-plate CV from sample control wells |
| `03_DA_v2.R` | Differential abundance per target, sensitivity refits, QC annotation |
| `04_QC_v2.R` | Plate ICC, PCA covariate screen, zero report, QC figures |
| `05_combine_da_results.R` | Stack per-cohort results into one table per plate model |
| `06_meta_v2.R` | Meta-analysis across cohorts |
| `06b_interaction_v2.R` | Phenotype × sex and phenotype × age interactions, meta-analysed |

## Input data

- **NPQ:** NULISA Protein Quantification values as delivered by Alamar (internal-control and inter-plate-control normalised, log2 scale). No further normalisation is applied.
- **GP2 release metadata:** phenotype, age, sex and ancestry, joined by donor ID. Donors not yet in the release take these fields from the cohort manifest; the source is recorded per sample.
- **Cohort manifests:** map site sample identifiers to GP2 donor IDs.

Input data are controlled-access and are not included in this repository. Results are written to `PROJECT_DIR` (set in the config).

## Model

```
NPQ ~ phenotype + age + sex (+ plate)
```

Fitted per target with `limma::lmFit` without empirical Bayes moderation, so ordinary standard errors are available for meta-analysis. Benjamini–Hochberg FDR < 0.05.

**Plate:** included as a covariate where both cases and controls are represented across plates. Left out where the within-plate comparison would rest on a small subset of cases. Every cohort is also run without plate (`plate_none`). 

Two sensitivity refits accompany each result: one without samples carrying an assay QC warning, one without PCA outliers (> 4 SD from the centroid of PCs 1–5).

## Meta-analysis

`05_combine_da_results.R` writes `DA_all_cohorts_combined_<model>.csv` for `plate_adj` (per-cohort plate rule) and `plate_none`. `06_meta_v2.R` reads the one named by `META_INPUT` and the cohorts in `META_COHORTS`, and fits per target with `metafor::rma`:

- Random effects (REML), z-test: primary
- Fixed effect: alongside
- BH FDR; I² and τ for heterogeneity; leave-one-cohort-out direction check

## Quality annotations

- **Detectability:** share of samples above the target's limit of detection. Targets below 50% (Alamar's threshold) are marked not `reportable`; high-abundance targets, which have no LOD, are exempt.
- **CV:** intra- and inter-plate, from sample control wells, on unlogged NPQ.
- **Zero imbalance:** difference in exact-zero rate between cases and controls.
- **Vendor flags** (Alamar Neuro220 data sheet): cross-reactivity, low vendor detectability, high vendor CV. Hemolysis-sensitive targets are flagged from Alamar's Inflammation 250 validation data (assumed to carry over to Neuro220; not vendor-confirmed).

## Environment

R 4.5. See `environment.yml`. Key packages: SummarizedExperiment, limma, metafor, lme4, sva, tidyverse.

## Status

Under active development as more cohorts are processed. Results are preliminary.
