# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# TO DO

FIND BIBLIOGRAPHY ON LIPIDOMICS, to undestand which normalization should I use, and which tests should I make ? 
Faire un rapport plus concis qui explique tous les resultats (enlever enrichment analysis etc, peut etre le gars pour une analyse annexe mais c'est pas l'important)

PLASMA/PBMC: 
- Regarder pour les lipides significatifs pour PLASMA a D0, et comparer boxplot avec resultat PBMC ( mettre les box plot cote a cote) ( ET VICE VERSA )
- faire la meme etude avec normalisation PQN pour voir les differences.
- Refaire analyse D60: Normaliser D0 par D60 pour chaque individu, et apres tu compare les lipides mild vs severe. (pour essayer d'enlever les biais individuels) . Est-ce que on enleve  la normalisation PQN, et on fait direct sur les valeurs brut? 
LONGITUDINAL:
- Regarder evolution au cours du temps, pour tous les lipides significatifs de l'etude PLASMA (que fait TIME X Severity acutellement ?) ; rajouter des stats, genre significativites.
- Regarder pour les lipides significatifs pour PBMC, et voir leur evolution dans PLASMA.

## Project overview

MISSE Global Lipidomics — a lipid-profiling study of dengue patients (**healthy** / **mild** / **severe**) across timepoints D0 (enrollment), D3, D10, D60 (convalescence). **The purpose of the study is to identify potential lipid biomarkers linked with dengue severity.** Sample naming convention: `PATIENT-ID` + `DXX`, e.g. `KT-313D0`, `BS-082D10`.

**Two distinct sample types exist in this repo`Dengue_lipidomics`:**
- **Plasma** — the mature, multi-timepoint dataset under `data/lipidsig_datasets/`. Used by essentially all existing analysis scripts (`01_analysis_three_groups*.R`, `02/03/04_posthoc_*.R`, `code/05_longitudinal_analysis.R`).
- **PBMC** — a newer, currently D0-only dataset under `data/PBMC/`, built by `code/PBMC_lipidomics_pretreatment.R` from `data/raw_data/260204_PBMC_DMI_Global_Results.xlsx`. Not yet documented in `README.md`.

Do not assume sample type from a file/directory name alone (e.g. `PBMC_01_analysis_three_groups.R` actually reads the *plasma* dataset) — check which folder under `data/` a script actually reads from.

**Clinical/study-design context:**
- **All `mild`/`severe` (sick) patients were hospitalized** — every non-healthy sample comes from a hospital admission. `healthy` samples are blood donor controls, not hospital patients.
- **`Day of Fever`** (a small integer, e.g. 0-7; not to be confused with the neighboring `Date of fever` column, which is a calendar date) lives in `data/raw_data/Final_Sample Shipment_ Montpellier_patient list_Dorothee_08082024 (2).xlsx` — default/first sheet ("Longitudinal samples"), header on row 2 (`openxlsx::read.xlsx(..., startRow = 2)`), columns `Patient.code`/`Day.of.Fever`. It records how many days the patient had already had fever *before* hospitalization/sampling. This means the D0 timepoint is enrollment/admission day, **not** a fixed point relative to infection or symptom onset — two patients both sampled at "D0" can be at different points in their actual illness course (severe patients in this cohort average ~1 day later presentation than mild: 3.24 vs 2.32 days).
  - A ready-to-use per-patient extract lives at `data/sick_patients_day_of_fever.tsv` (43 patients, matches the full longitudinal cohort). It's produced by `code/PBMC_lipidomics_pretreatment.R` (search `sick_patients_fever`) — **that script previously had a bug** selecting `Date.of.fever` (the calendar date) instead of `Day.of.Fever` and mislabeling it; fixed 2026-07-13, regenerate via that script section if the source xlsx changes.
  - Used as a covariate in `05_longitudinal_analysis.R`'s LMM (`value ~ timepoint * severity + day_of_fever + (1 | patient_id)`, present in every nested model so it never changes what an LRT isolates). In practice this barely moved results (same 23 interaction-significant lipids, nearly identical p-values) — expected, since a per-patient-constant covariate is largely redundant with what the `(1 | patient_id)` random intercept already absorbs. A stronger test of "is D0 comparable across patients" would need to realign each patient's timeline to days-since-fever-onset rather than just adding this as a fixed-effect covariate.

`README.md` has the full repository structure, step-by-step workflow, and key parameter tables for the plasma cross-sectional pipeline — read it for file-level detail not repeated here.

## Environment & commands

- R 4.5.0, invoked as `C:\Program Files\R\R-4.5.0\bin\Rscript.exe` on this Windows machine (not on PATH in git-bash; use the full path, or PowerShell).
- renv-managed. **The active lockfile is `renv.lock` at the repo root**, loaded automatically via `.Rprofile` (`source("renv/activate.R")`). The `LipidSigR/renv/` subfolder is a leftover macOS/ARM copy from the original author's machine (contains `R-4.3/aarch64-apple-darwin20` paths) and is *not* the live environment — ignore it, despite `README.md` §3 pointing there.
- Restore packages: `Rscript -e "renv::restore()"` from the repo root. A `renv::status()` "project is out-of-sync" message printed on every `Rscript` invocation is a known, harmless artifact here — scripts still run correctly.
- Run an analysis script end-to-end: `Rscript code/<script>.R`.
- Render a markdown/Rmd report to HTML: `Rscript -e "rmarkdown::render('path/to/Report.Rmd', output_format='html_document')"`. Match the existing house style — YAML `html_document` with `theme: flatly`, `toc_float: true`, `number_sections: true` (see `HOWTOGENERATEREPORT.md` and any `analysis/**/*.Rmd`/`.md` report). With `number_sections: true`, do not also hand-number markdown headers (`## 1. Foo`) — pandoc numbers them automatically and the two schemes double up.
- Python preprocessing notebooks (`code/*.ipynb`) need `pandas numpy matplotlib seaborn scipy openpyxl`, run via Jupyter — see README §3 for the pip/conda install line.

## Architecture

- **Two-stage pipeline.** Python notebooks turn raw platform Excel exports (`data/raw_data/`) into per-timepoint TSVs (`data/pretreatment/`), which are merged into LipidSigR-ready `(abundance, group_information)` TSV pairs (`data/lipidsig_datasets/`, or `data/PBMC/` for the PBMC arm). R scripts only ever consume the latter, never the raw Excel files directly.
- **LipidSigR + rgoslin do the heavy lifting.** Standard call sequence: `rgoslin::parseLipidNames()` (drop `Grammar == "NOT_PARSEABLE"` rows) → `as_summarized_experiment()` → `data_process()` → profiling/DE/enrichment functions. Standard `data_process()` parameters used throughout: `exclude_missing_pct=70`, `replace_na_method="min"` (×0.5 ref), `normalization="Percentage"`, `transform="log10"`. On that scale a value is `log10(% of total lipidome in that sample)` — negative numbers are expected for any lipid under 1% of the total, not a bug.
- **Cross-sectional pipeline** (`01_analysis_three_groups*.R` → `02/03/04_posthoc_*.R`): one timepoint at a time. 3-group ANOVA (`deSp_multiGroup`, `se_type="de_multiple"`), followed by 3 pairwise 2-group t-tests (`deSp_twoGroup`, `se_type="de_two"`, ref group `healthy`/`mild`). Each writes a parallel output tree under `analysis/<Comparison>/<Timepoint>/` (`01.Profiling/`, `02.DiffExp/`, `03.Enrichment/`).
- **Longitudinal pipeline** (`code/05_longitudinal_analysis.R`): tracks the same 43 mild/severe patients across all 4 timepoints. Deliberately does **not** use `deSp_multiGroup`'s ANOVA/Kruskal-Wallis for the time comparison — treating 4 repeated measurements per patient as independent samples is pseudoreplication. Instead fits a per-lipid linear mixed model (`lme4::lmer`, `value ~ timepoint * severity + day_of_fever + (1 | patient_id)`, LRT-based p-values via nested-model comparison, BH-FDR correction), reusing LipidSigR only for annotation/QC/normalization/PCA (all per-sample operations that don't assume independence). The same LMM approach is reused for class- and category-level aggregation via rgoslin's own `Lipid.Maps.Main.Class` / `Lipid.Maps.Category` columns (species values are summed in linear space, not log space, before re-aggregating).
- `data/lipidsig_datasets/Longitudinal_matched/Lipid_abundance_data_{D0,D03,D10,D60}.tsv` (renamed 2026-07-24 from `Deprecated/`, which was misleading — it's the only dataset where the same patients are matched by raw abundance across all 4 timepoints in consistent column order) is what the longitudinal pipeline reads. Other "current" folders (`healthy_vs_sick_patients/`, `selected_patients/`) each cover only 2-3 of the 4 timepoints. The other ~16 files that used to sit alongside these in `Deprecated/` (Log2FC variants, `*_severe.tsv` pairwise files, an internally-inconsistent `D0_vs_D60_severe` draft) were unreferenced by any script and were deleted rather than moved.
- **Known LipidSigR bug:** `as_summarized_experiment()` throws `invalid rownames length` if the input has both constant (all-zero) features *and* constant samples at once — its internal constant-value filter desyncs the annotation table across the two axes. Pre-filter constant rows/columns yourself before calling it (see `05_longitudinal_analysis.R`, Data Loading section, for the pattern).
- **Hardcoded paths:** most existing R scripts hardcode the original author's macOS path (`/Users/loictalignani/research/project/lipidomics`) for `setwd()` and the output `base_dir` — must be edited per-machine (README §6 has a sed one-liner). `05_longitudinal_analysis.R` instead derives everything from `getwd()`/relative paths — prefer that pattern in new scripts rather than another hardcoded absolute path.
- `code/CHIKV_lipidomics_analysis.py` is an unrelated in-vitro sub-study (CHIKV drug treatment in Vero/C636 cell lines) — not part of the dengue patient cohort, safe to ignore for biomarker-discovery work.
