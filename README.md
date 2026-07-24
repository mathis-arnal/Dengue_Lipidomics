---
editor_options: 
  markdown: 
    wrap: 72
---

# MISSE Global Lipidomics — Analysis Pipeline

Lipid profiling study of DENGUE patients across multiple time points
(D0, D3, D10, D60), comparing three clinical groups: **healthy**,
**mild**, and **severe**.

**Contact:** Loïc Talignani —
[loic.talignani\@ird.fr](mailto:loic.talignani@ird.fr){.email}

------------------------------------------------------------------------

## Table of Contents

1.  [Study Design](#1-study-design)

2.  [Repository Structure](#2-repository-structure)

3.  [Software Requirements](#3-software-requirements)

4.  [Analysis Workflow — Step by
    Step](#4-analysis-workflow--step-by-step)

    -   [Step 1 — Raw data to TSV
        (Python)](#step-1--raw-data-to-tsv-python)
    -   [Step 2 — Merge & prepare datasets
        (Python)](#step-2--merge--prepare-datasets-python)
    -   [Step 3 — Three-group analysis
        (R)](#step-3--three-group-analysis-r)
    -   [Step 4 — Post-hoc pairwise comparisons
        (R)](#step-4--post-hoc-pairwise-comparisons-r)

5.  [Output Structure](#5-output-structure)

6.  [Adapting to a New Machine](#6-adapting-to-a-new-machine)

7.  [Key Parameters](#7-key-parameters)

------------------------------------------------------------------------

## 1. Study Design

| Group     | Description                     |
|-----------|---------------------------------|
| `healthy` | Healthy blood donors (controls) |
| `mild`    | DENGUE patients — mild form     |
| `severe`  | DENGUE patients — severe form   |

**Time points:** D0 (enrollment/baseline), D3 (day 3), D10 (day 10), D60
(day 60 — convalescence).

**Sample naming convention:** `PATIENT-ID` + `DXX` e.g. `KT-313D0`,
`BS-082D10`.

The primary analysis compares the three groups at D0 and D3, using an
ANOVA-like approach followed by pairwise post-hoc tests (Healthy vs
Mild, Healthy vs Severe, Mild vs Severe).

------------------------------------------------------------------------

## 2. Repository Structure

```         
lipidomics/
│
├── code/                                   # All analysis scripts
│   │
│   │   # Python preprocessing notebooks
│   ├── lipidomics.ipynb                    # Initial preprocessing: Excel → per-timepoint TSV
│   ├── healthy\\\_vs\\\_sick\\\_patients\\\_pretreatment.ipynb   # Merge healthy + sick patients (D0)
│   ├── healthy\\\_vs\\\_sick\\\_patients\\\_pretreatment\\\_D3.ipynb # Same for D3
│   ├── Log2FC\\\_calculation.ipynb            # Log2 fold-change computation (individual lipids)
│   ├── Log2FC\\\_calculation\\\_sum\\\_species.ipynb # Log2FC on summed lipid species
│   ├── Normalization\\\_by\\\_division.ipynb     # Ratio-based normalization
│   ├── Normalization\\\_by\\\_division\\\_v2.ipynb  # Updated normalization
│   ├── Oxylipins\\\_pretreatment.ipynb        # Oxylipin-specific preprocessing
│   ├── Oxylipins\\\_pretreatment\\\_v2.ipynb     # Updated oxylipin preprocessing
│   ├── chi2\\\_analysis\\\_ratios.ipynb          # Chi-squared test on lipid ratios
│   ├── chi2\\\_analysis\\\_sum\\\_species.ipynb     # Chi-squared test on lipid species sums
│   │
│   │   # R — Main LipidSigR analysis scripts (three groups + post-hoc)
│   ├── 01\\\_analysis\\\_three\\\_groups.R          # ANOVA / 3-group analysis at D0
│   ├── 01\\\_analysis\\\_three\\\_groups\\\_D3.R       # Same analysis at D3
│   ├── 02\\\_posthoc\\\_healthy\\\_vs\\\_mild.R        # Post-hoc: Healthy vs Mild (D0)
│   ├── 02\\\_posthoc\\\_healthy\\\_vs\\\_mild\\\_D3.R     # Post-hoc: Healthy vs Mild (D3)
│   ├── 03\\\_posthoc\\\_healthy\\\_vs\\\_severe.R      # Post-hoc: Healthy vs Severe (D0)
│   ├── 03\\\_posthoc\\\_healthy\\\_vs\\\_severe\\\_D3.R   # Post-hoc: Healthy vs Severe (D3)
│   ├── 04\\\_posthoc\\\_mild\\\_vs\\\_severe.R         # Post-hoc: Mild vs Severe (D0)
│   ├── 04\\\_posthoc\\\_mild\\\_vs\\\_severe\\\_D3.R      # Post-hoc: Mild vs Severe (D3)
│   ├── profiling\\\_posthoc\\\_comparisons.R     # Profiling plots across comparisons
│   ├── generate\\\_report.R                   # Report generation helper
│   │
│   │   # R — Legacy / complementary analyses
│   ├── LipidSigR\\\_analysis\\\_selected\\\_patients.R    # 2-group analysis (selected patients, D0)
│   ├── LipidSigR\\\_analysis\\\_selected\\\_patients\\\_D3.R # Same at D3
│   ├── LipidSigR\\\_analysis\\\_selected\\\_patients\\\_D10.R # Same at D10
│   ├── LipidSigR\\\_analysis\\\_Healthy\\\_vs\\\_sick\\\_patients.R # Healthy vs sick (2-group)
│   ├── LipidSigR\\\_D0\\\_D60.R                 # Paired analysis D0 vs D60 (severe patients only)
│   │
│   │   # Python — CHIKV in-vitro analysis
│   └── CHIKV\\\_lipidomics\\\_analysis.py        # In-vitro CHIKV drug treatment analysis (Vero/C636)
│
├── data/
│   ├── raw\\\_data/                           # Original files from the lipidomics platform
│   │   ├── MISSE\\\_Global\\\_Lipidomic\\\_results.xlsx          # Main raw data (all patients, all timepoints)
│   │   ├── 241008GMI\\\_oxylipins\\\_results.xlsx             # Oxylipin raw data (batch 1)
│   │   ├── 241009DMI\\\_oxylipins\\\_results.xlsx             # Oxylipin raw data (batch 2)
│   │   ├── MISSE\\\_200800835\\\_2025\\\_May\\\_donneurs\\\_sains/      # Healthy donor raw data (2025 batch)
│   │   │   └── 250528DMI\\\_Global\\\_Lipidomic\\\_results.xlsx
│   │   ├── Demographic\\\_data.xlsx                        # Patient demographic information
│   │   └── Patient\\\_data.xlsx                            # Patient clinical data
│   │
│   ├── pretreatment/                       # Intermediate files produced by Python notebooks
│   │   ├── MISSE\\\_Global\\\_Lipidomics\\\_results\\\_all.tsv       # All patients, all timepoints (TSV)
│   │   ├── MISSE\\\_Global\\\_Lipidomics\\\_results\\\_all.csv       # Same in CSV
│   │   ├── MISSE\\\_Global\\\_Lipidomics\\\_results\\\_all\\\_with\\\_log2FC.tsv  # With log2FC columns
│   │   └── healthy\\\_vs\\\_sick\\\_patients/                    # Merged healthy + patient datasets
│   │       ├── Healthy\\\_patients.tsv                     # Healthy donor abundance data
│   │       ├── MISSE\\\_Global\\\_Lipidomics\\\_results\\\_selected\\\_patients\\\_D0.tsv
│   │       ├── MISSE\\\_Global\\\_Lipidomics\\\_results\\\_selected\\\_patients\\\_D3.tsv
│   │       ├── MISSE\\\_Global\\\_Lipidomics\\\_results\\\_selected\\\_patients\\\_D10.tsv
│   │       ├── healthy\\\_sick\\\_lipidomics.tsv              # Merged abundance (D0)
│   │       └── healthy\\\_sick\\\_lipidomics\\\_D3.tsv           # Merged abundance (D3)
│   │
│   └── lipidsig\\\_datasets/                  # Final formatted datasets — input to R scripts
│       ├── healthy\\\_vs\\\_sick\\\_patients/        # PRIMARY DATASETS (used by 01-04 scripts)
│       │   ├── healthy\\\_sick\\\_lipidomics.tsv             # Abundance: all 3 groups, D0
│       │   ├── healthy\\\_sick\\\_lipidomics\\\_D3.tsv          # Abundance: all 3 groups, D3
│       │   ├── healthy\\\_MILD\\\_lipidomics.tsv             # Abundance: healthy + mild, D0
│       │   ├── healthy\\\_MILD\\\_lipidomics\\\_D3.tsv          # Abundance: healthy + mild, D3
│       │   ├── healthy\\\_SEVERE\\\_lipidomics.tsv           # Abundance: healthy + severe, D0
│       │   ├── healthy\\\_SEVERE\\\_lipidomics\\\_D3.tsv        # Abundance: healthy + severe, D3
│       │   ├── MILD\\\_SEVERE\\\_lipidomics.tsv              # Abundance: mild + severe, D0
│       │   ├── MILD\\\_SEVERE\\\_lipidomics\\\_D3.tsv           # Abundance: mild + severe, D3
│       │   ├── group\\\_information\\\_table\\\_healthy\\\_vs\\\_sick\\\_patients\\\_D0.tsv  # Groups for 3-group
│       │   ├── group\\\_information\\\_table\\\_healthy\\\_vs\\\_sick\\\_patients\\\_D3.tsv
│       │   ├── group\\\_information\\\_table\\\_healthy\\\_vs\\\_MILD\\\_D0.tsv
│       │   ├── group\\\_information\\\_table\\\_healthy\\\_vs\\\_MILD\\\_D3.tsv
│       │   ├── group\\\_information\\\_table\\\_healthy\\\_vs\\\_SEVERE\\\_D0.tsv
│       │   ├── group\\\_information\\\_table\\\_healthy\\\_vs\\\_SEVERE\\\_D3.tsv
│       │   ├── group\\\_information\\\_table\\\_MILD\\\_vs\\\_SEVERE\\\_patients\\\_D0.tsv
│       │   └── group\\\_information\\\_table\\\_MILD\\\_vs\\\_SEVERE\\\_patients\\\_D3.tsv
│       ├── selected\\\_patients/               # Subset of selected patients only (no healthy)
│       ├── Oxylipines/                      # Oxylipin-specific datasets
│       ├── Ratios/                          # Lipid ratio datasets
│       └── Longitudinal\\\_matched/            # Same 43 mild/severe patients, raw abundances, all 4 timepoints (used by 05/07/08 scripts)
│
├── analysis/                               # OUTPUT — generated figures and tables
│   ├── Three\\\_groups/
│   │   ├── D0/                             # Results for D0 three-group analysis
│   │   └── D3/                             # Results for D3 three-group analysis
│   ├── PostHoc\\\_Healthy\\\_vs\\\_Mild/
│   │   ├── D0/
│   │   └── D3/
│   ├── PostHoc\\\_Healthy\\\_vs\\\_Severe/
│   │   ├── D0/
│   │   └── D3/
│   ├── PostHoc\\\_Mild\\\_vs\\\_Severe/
│   │   ├── D0/
│   │   └── D3/
│   ├── Healthy\\\_vs\\\_sick\\\_patients/           # Legacy 2-group outputs
│   ├── Selected\\\_patients/                  # Legacy selected-patient outputs
│   ├── CHIKV/                              # In-vitro CHIKV analysis outputs
│   └── deprecated/                         # Old analyses
│
├── LipidSigR/                              # Local R environment managed by renv
│   └── renv.lock                           # Package lockfile — do not modify manually
│
├── PPTs/                                   # Presentation slides and reports
├── biblio/                                 # Bibliography
└── CLAUDE.md                               # Instructions for Claude Code assistant
```

------------------------------------------------------------------------

## 3. Software Requirements

### Python (preprocessing)

-   Python 3.8+
-   `pandas`, `numpy`, `matplotlib`, `seaborn`, `scipy`, `openpyxl`

Install with:

``` bash
pip install pandas numpy matplotlib seaborn scipy openpyxl
```

Or with conda:

``` bash
conda install pandas numpy matplotlib seaborn scipy openpyxl
```

Notebooks are run in **Jupyter Lab** or **Jupyter Notebook**.

### R (analysis)

Core packages:

-   **LipidSigR** — from GitHub: `BioinfOMICS/LipidSigR`
-   **rgoslin** — Bioconductor, for lipid name standardization
-   **SummarizedExperiment** — Bioconductor
-   **tidyverse** (dplyr, readr, ggplot2)
-   **ggplot2**, **ggsignif** — for post-hoc visualizations

Full installation:

``` r
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("rgoslin", "SummarizedExperiment", "fgsea",
                       "gatom", "mixOmics", "S4Vectors", "BiocGenerics"))

# LipidSigR from GitHub
install.packages("devtools")
devtools::install\\\_github("BioinfOMICS/LipidSigR",
                         build\\\_vignettes = TRUE,
                         dependencies = TRUE)

# CRAN packages
install.packages(c("tidyverse", "ggplot2", "ggsignif", "plotly",
                   "factoextra", "ggthemes", "ggforce", "Hmisc",
                   "heatmaply", "Rtsne", "uwot", "rsample",
                   "ranger", "caret", "yardstick"))
```

**Recommended:** use the `renv` lockfile in `LipidSigR/` for a
reproducible environment:

``` r
setwd("LipidSigR/")
renv::restore()   # Installs exact package versions from lockfile
```

------------------------------------------------------------------------

## 4. Analysis Workflow — Step by Step

### Step 1 — Raw data to TSV (Python)

**Script:** `code/lipidomics.ipynb`\
**Input:**
`data/raw\\\_data/MISSE\\\_Global\\\_Lipidomic\\\_results.xlsx`\
**Output:**
`data/pretreatment/MISSE\\\_Global\\\_Lipidomics\\\_results\\\_all.tsv`
(and per-timepoint TSV files)

The raw Excel file from the lipidomics platform contains all patients
across all time points in a single sheet. Values use **commas as decimal
separators** (European format). This notebook:

1.  Reads the Excel file with `pandas`
2.  Replaces commas with periods in all numeric columns
3.  Splits data by time point (D0, D3, D10, D60) using column name
    suffix filtering
4.  Saves each time point as a separate TSV

**TSV format produced:**

| feature | PATIENT-ID1 | PATIENT-ID2 | …   |
|---------|-------------|-------------|-----|
| DG 32:0 | 0.75        | 1.84        | …   |
| TG 54:3 | 3.21        | 0.00        | …   |

-   First column is always `feature` (lipid name)
-   Subsequent columns are sample IDs with time point suffix (e.g.
    `KT-313D0`)
-   Missing values are filled with `0`

------------------------------------------------------------------------

### Step 2 — Merge & prepare datasets (Python)

**Scripts:**

-   `code/healthy\\\_vs\\\_sick\\\_patients\\\_pretreatment.ipynb` (for
    D0)
-   `code/healthy\\\_vs\\\_sick\\\_patients\\\_pretreatment\\\_D3.ipynb`
    (for D3)

**Purpose:** Merge the healthy donor dataset with the sick patient
dataset, and produce the two required LipidSigR input files: **abundance
table** and **group information table**.

#### 2a. Abundance table

1.  Load healthy donors:
    `data/pretreatment/healthy\\\_vs\\\_sick\\\_patients/Healthy\\\_patients.tsv`
2.  Load sick patients:
    `data/pretreatment/healthy\\\_vs\\\_sick\\\_patients/MISSE\\\_Global\\\_Lipidomics\\\_results\\\_selected\\\_patients\\\_D0.tsv`
3.  Clean lipid names (replace commas in feature names with semicolons)
4.  Merge on `feature` column with `inner join` (keeps only lipids
    present in both datasets)
5.  Save to:
    `data/lipidsig\\\_datasets/healthy\\\_vs\\\_sick\\\_patients/healthy\\\_sick\\\_lipidomics.tsv`

The notebook also produces subset files for each post-hoc pair:

-   `healthy\\\_MILD\\\_lipidomics.tsv` — healthy + mild patients only
-   `healthy\\\_SEVERE\\\_lipidomics.tsv` — healthy + severe patients
    only
-   `MILD\\\_SEVERE\\\_lipidomics.tsv` — mild + severe patients only

#### 2b. Group information table

This file tells LipidSigR which group each sample belongs to.

Required format (TSV, 4 columns):

| sample_name      | label_name       | group   | pair |
|------------------|------------------|---------|------|
| KT-313D0         | KT-313D0         | mild    | NA   |
| KT-193D0         | KT-193D0         | severe  | NA   |
| JV-015-RAN-IC-M2 | JV-015-RAN-IC-M2 | healthy | NA   |

-   `sample\\\_name`: must match column headers in the abundance table
    exactly
-   `group`: one of `healthy`, `mild`, `severe`
-   `pair`: `NA` for unpaired designs

Files are saved to:
`data/lipidsig\\\_datasets/healthy\\\_vs\\\_sick\\\_patients/`

------------------------------------------------------------------------

### Step 3 — Three-group analysis (R)

**Scripts:**

-   `code/01\\\_analysis\\\_three\\\_groups.R` — D0
-   `code/01\\\_analysis\\\_three\\\_groups\\\_D3.R` — D3

**Input files:**

```         
data/lipidsig\\\_datasets/healthy\\\_vs\\\_sick\\\_patients/healthy\\\_sick\\\_lipidomics.tsv
data/lipidsig\\\_datasets/healthy\\\_vs\\\_sick\\\_patients/group\\\_information\\\_table\\\_healthy\\\_vs\\\_sick\\\_patients\\\_D0.tsv
```

**Output directory:** `analysis/Three\\\_groups/D0/` (or `D3/`)

#### What this script does

1.  **Data loading** — reads abundance + group info with `read\\\_tsv()`

2.  **Lipid name parsing** — standardizes nomenclature with
    `rgoslin::parseLipidNames()`

    -   Lipids not recognized (`Grammar == "NOT\\\_PARSEABLE"`) are
        filtered out

3.  **SummarizedExperiment object** — built with
    `as\\\_summarized\\\_experiment(..., se\\\_type = "de\\\_multiple")`

4.  **Data processing** via `data\\\_process()`:

    -   Missing value filter: keep lipids present in ≥ 70% of samples
    -   NA replacement: minimum value × 0.5
    -   Normalization: `"Percentage"`
    -   Transformation: `log10`

5.  **Profiling** (saved to `01.Profiling/`):

    -   Data quality box plots and density plots (before/after
        processing)
    -   Cross-sample variability plots
    -   Dimensionality reduction: **PCA**, **t-SNE**, **UMAP**
    -   Correlation heatmaps (by sample, by lipid category, by lipid
        class)
    -   Lipid characteristics plots

6.  **Differential expression** — ANOVA across the three groups
    (`deSp\\\_multiGroup()`), saved to `02.DiffExp/`

7.  **Enrichment analysis** — LSEA and ORA, saved to `03.Enrichment/`

#### How to run

``` r
# In RStudio or R terminal — update the path first (see Section 6)
source("code/01\\\_analysis\\\_three\\\_groups.R")
```

------------------------------------------------------------------------

### Step 4 — Post-hoc pairwise comparisons (R)

After the three-group ANOVA, three pairwise comparisons are performed.
Each uses a **two-group t-test** (or Wilcoxon) and is run independently.

| Script | Comparison | Input abundance file |
|----|----|----|
| `02\\\_posthoc\\\_healthy\\\_vs\\\_mild.R` | Healthy vs Mild | `healthy\\\_MILD\\\_lipidomics.tsv` |
| `03\\\_posthoc\\\_healthy\\\_vs\\\_severe.R` | Healthy vs Severe | `healthy\\\_SEVERE\\\_lipidomics.tsv` |
| `04\\\_posthoc\\\_mild\\\_vs\\\_severe.R` | Mild vs Severe | `MILD\\\_SEVERE\\\_lipidomics.tsv` |

D3 variants follow the same logic (`\\\_D3.R` suffix).

#### What each post-hoc script does

1.  **Data loading & processing** — identical pipeline to Step 3 but
    with `se\\\_type = "de\\\_two"`

2.  **Profiling** — PCA, t-SNE, UMAP, correlation heatmaps (saved to
    `01.Profiling/`)

3.  **Differential expression** — `deSp\\\_twoGroup()` with:

    -   Test: `t-test` (or `wilcox`)
    -   p-value cutoff: `0.05`
    -   Fold-change cutoff: `1.5`
    -   Reference group: `healthy` (or `mild` for mild vs severe)

4.  **Visualizations** saved to `02.DiffExp/02.Visualizations/`:

    -   Volcano plot
    -   MA plot
    -   Heatmap of significant lipids

5.  **Individual lipid boxplots** for each significant lipid

6.  **Enrichment analysis** — LSEA + ORA (saved to `03.Enrichment/`)

#### How to run (all post-hoc comparisons at D0)

``` r
source("code/02\\\_posthoc\\\_healthy\\\_vs\\\_mild.R")
source("code/03\\\_posthoc\\\_healthy\\\_vs\\\_severe.R")
source("code/04\\\_posthoc\\\_mild\\\_vs\\\_severe.R")
```

**Recommended order:** always run `01\\\_analysis\\\_three\\\_groups.R`
before the post-hoc scripts, as the ANOVA identifies which lipids to
follow up on.

------------------------------------------------------------------------

## 5. Output Structure

Each analysis directory follows the same sub-structure:

```         
analysis/Three\\\_groups/D0/
├── 01.Profiling/
│   ├── 00.Data\\\_quality/
│   │   ├── BoxPlot\\\_before\\\_process.png
│   │   ├── BoxPlot\\\_after\\\_process.png
│   │   ├── DensityPlot\\\_before\\\_process.png
│   │   └── DensityPlot\\\_after\\\_process.png
│   ├── 01.Cross-sample\\\_variability/
│   │   ├── Expressed\\\_lipid\\\_numbers.png
│   │   ├── Histogram\\\_lipid\\\_amount.png
│   │   └── Lipid\\\_abundance\\\_distribution.png
│   ├── 02.Dimensionality\\\_reduction/
│   │   ├── PCA/    → PCA\\\_plot.png, loadings, scree plot
│   │   ├── t-SNE/  → tSNE\\\_plot.png
│   │   └── UMAP/   → UMAP\\\_plot.png
│   ├── 03.Correlation\\\_Heatmap/
│   │   ├── by\\\_samples/
│   │   ├── by\\\_category/
│   │   └── by\\\_class/
│   └── 04.Lipid\\\_characteristics/
├── 02.DiffExp/
│   ├── 01.ANOVA/               # (or 01.t-test\\\_results/ for post-hoc)
│   ├── 02.Visualizations/
│   │   ├── Volcano/
│   │   ├── MA\\\_plot/
│   │   └── Heatmap/
│   └── 03.Individual\\\_lipid\\\_boxplots/
└── 03.Enrichment/
    ├── LSEA/
    └── ORA/
```

------------------------------------------------------------------------

## 6. Adapting to a New Machine

All R scripts contain hardcoded absolute paths that must be updated
before running:

``` r
# Find and replace in all .R files:
setwd("/Users/loictalignani/research/project/lipidomics")
base\\\_dir <- "/Users/loictalignani/research/project/lipidomics/analysis/..."
```

**Replace** `/Users/loictalignani/research/project/lipidomics` with your
own project path.

On Linux/macOS, a quick way to update all scripts at once:

``` bash
# Preview changes first
grep -r "loictalignani" code/\\\*.R

# Replace (macOS sed syntax)
sed -i '' 's|/Users/loictalignani/research/project/lipidomics|/YOUR/NEW/PATH|g' code/\\\*.R

# Linux sed syntax
sed -i 's|/Users/loictalignani/research/project/lipidomics|/YOUR/NEW/PATH|g' code/\\\*.R
```

------------------------------------------------------------------------

## 7. Key Parameters

### Data processing (`data\\\_process`)

| Parameter | Value used | Description |
|----|----|----|
| `exclude\\\_missing\\\_pct` | `70` | Keep lipids present in ≥ 70% of samples |
| `replace\\\_na\\\_method` | `"min"` | Replace missing values with minimum |
| `replace\\\_na\\\_method\\\_ref` | `0.5` | Multiply minimum by 0.5 |
| `normalization` | `"Percentage"` | Normalize to percentage of total |
| `transform` | `"log10"` | Log10 transformation |

### Differential expression

| Script type | Function               | `se\\\_type`       | Test   |
|-------------|------------------------|--------------------|--------|
| Three-group | `deSp\\\_multiGroup()` | `"de\\\_multiple"` | ANOVA  |
| Post-hoc    | `deSp\\\_twoGroup()`   | `"de\\\_two"`      | t-test |

### Lipid name parsing (rgoslin)

Grammar used: `"FattyAcids"` (default) — covers most standard lipidomics
nomenclature. Alternative grammars: `"Shorthand2020"`, `"Goslin"`,
`"LipidMaps"`, `"SwissLipids"`, `"HMDB"`.

Lipids with `Grammar == "NOT\\\_PARSEABLE"` are excluded from all
downstream analyses. Typical recognition rate: \~85–95% of detected
features.
