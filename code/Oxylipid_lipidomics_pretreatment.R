# Oxylipid lipidomics — pretreatment
#
# Builds per-timepoint (D0/D3/D10/D60) abundance + group-information TSVs for the
# oxylipid panel, mirroring PBMC_lipidomics_pretreatment.R's pattern.
#
# Source: data/raw_data/241008GMI_oxylipins_results.xlsx, sheet "Results" — the
# platform-reported, already-deduplicated, already-normalized (pg/uL of fluid)
# results table: 47 oxylipins with a single value per compound per sample, plus an
# authoritative parent-PUFA "Precursor" label per compound (column 1, filled down).
#
# This supersedes an earlier version of this script that read
# data/raw_data/oxylipines.xlsx (a pre-existing, undocumented 24-compound extract).
# That file turned out to omit roughly half of what the platform actually measured —
# missing essentially the entire specialized pro-resolving mediator class (resolvins
# RvE1/RvD1/RvD2/RvD3/RvD5, protectins PD1/PDx, maresin 7(S)-MaR1), several more
# prostaglandins/leukotrienes/lipoxins, and 3 of 4 EET regioisomers. Reading directly
# from the platform's own "Results" sheet recovers the full 47-compound panel.
#
# Two near-identical raw files exist (241008GMI_oxylipins_results.xlsx,
# 241009DMI_oxylipins_results.xlsx) despite the "batch 1"/"batch 2" naming in
# README.md — their "Raw Datas" sheets were confirmed byte-identical (16,830/16,830
# cells match), and their "Results" sheets differ only by 2 blank trailing columns in
# the DMI file. These are not two different sample batches; only one file is needed.
#
# There is no healthy arm for oxylipids: healthy blood donors are single-draw
# controls with no D60 sample, so this dataset is mild/severe only.

library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(readr)
library(zoo)

out_dir <- "data/lipidsig_datasets/Oxylipines_raw"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------
# 1. Load the Results sheet: 47 compounds x sample columns, + Precursor (family)
# ------------------------------------------------------------------------

raw_file <- "data/raw_data/241008GMI_oxylipins_results.xlsx"
results_grid <- openxlsx::read.xlsx(raw_file, sheet = "Results", colNames = FALSE)

# Row 4 holds sample IDs from column 3 onward; row 6 onward holds one row per
# compound (column 1 = precursor, filled down from group-start rows only;
# column 2 = compound name; columns 3+ = pg/uL-of-fluid values). Rows after the
# last compound are "Total oxylipins issued from <precursor> cascade" and a grand
# "Total" row -- platform-computed family sums, not individual features, dropped here.
sample_cols_idx <- 3:ncol(results_grid)
sample_ids <- trimws(as.character(unlist(results_grid[4, sample_cols_idx])))

precursor_filled <- zoo::na.locf(results_grid[[1]], na.rm = FALSE)
compound_name <- trimws(as.character(results_grid[[2]]))

is_compound_row <- !is.na(compound_name) &
  !grepl("^Total", compound_name) &
  compound_name != "pg/µL of fluid" &
  seq_len(nrow(results_grid)) >= 6

oxy_raw <- as.data.frame(results_grid[is_compound_row, sample_cols_idx])
colnames(oxy_raw) <- sample_ids
oxy_raw <- oxy_raw %>% mutate(across(everything(), as.numeric))
oxy_raw <- tibble(feature = compound_name[is_compound_row], oxy_raw)

precursor_by_feature <- tibble(
  feature = compound_name[is_compound_row],
  precursor = trimws(gsub("[\r\n]+", " ", precursor_filled[is_compound_row]))
)
write_tsv(precursor_by_feature, file.path(out_dir, "oxylipin_precursor_by_feature.tsv"))

cat(sprintf("Loaded %d compounds x %d samples from '%s' (Results sheet)\n",
            nrow(oxy_raw), ncol(oxy_raw) - 1, raw_file))

timepoints <- c(D0 = "D0", D03 = "D03", D10 = "D10", D60 = "D60")

split_by_timepoint <- function(tp_suffix) {
  sample_cols <- grep(paste0("^[A-Z]{2}-\\d+", tp_suffix, "$"), colnames(oxy_raw), value = TRUE)
  patient_ids <- sort(sub(paste0(tp_suffix, "$"), "", sample_cols))
  sample_cols <- paste0(patient_ids, tp_suffix)
  df <- oxy_raw %>% dplyr::select(feature, dplyr::all_of(sample_cols))
  df
}

abundance_by_tp <- lapply(timepoints, split_by_timepoint)

# ------------------------------------------------------------------------
# 2. Severity classification (mild/severe) from the master patient list
# ------------------------------------------------------------------------

sick_patients_df <- read.xlsx(
  "data/raw_data/Final_Sample Shipment_ Montpellier_patient list_Dorothee_08082024 (2).xlsx",
  startRow = 2)

sick_patients_status <- sick_patients_df %>%
  dplyr::select(Patient.code, `1997.Classification`) %>%
  tidyr::drop_na(`1997.Classification`) %>%
  mutate(status = ifelse(`1997.Classification` %in% c("DSS", "DHF"), "severe", "mild")) %>%
  dplyr::select(patient = Patient.code, status)

# ------------------------------------------------------------------------
# 3. Write abundance + group-information TSVs per timepoint
# ------------------------------------------------------------------------

for (tp_name in names(timepoints)) {
  tp_suffix <- timepoints[[tp_name]]
  abundance <- abundance_by_tp[[tp_name]]

  sample_cols <- setdiff(colnames(abundance), "feature")
  patient_ids <- sub(paste0(tp_suffix, "$"), "", sample_cols)

  group_info <- tibble(
    sample_name = sample_cols,
    patient = patient_ids
  ) %>%
    dplyr::left_join(sick_patients_status, by = "patient")

  unclassified <- group_info %>% dplyr::filter(is.na(status))
  if (nrow(unclassified) > 0) {
    cat(sprintf(
      "[%s] Dropping %d patient(s) with no severity classification in the master list: %s\n",
      tp_name, nrow(unclassified), paste(unclassified$patient, collapse = ", ")))
  }
  group_info <- group_info %>% tidyr::drop_na(status)

  abundance <- abundance %>% dplyr::select(feature, dplyr::all_of(group_info$sample_name))

  group_info_out <- tibble(
    sample_name = group_info$sample_name,
    label_name = group_info$sample_name,
    group = group_info$status,
    pair = NA_character_
  )

  abundance_file <- file.path(out_dir, sprintf("Oxylipid_lipid_abundance_data_%s.tsv", tp_name))
  group_file <- file.path(out_dir, sprintf("group_information_table_oxylipid_%s.tsv", tp_name))

  write_tsv(abundance, abundance_file)
  write_tsv(group_info_out, group_file)

  cat(sprintf(
    "[%s] %d features x %d patients written (%d mild, %d severe) -> %s / %s\n",
    tp_name, nrow(abundance), nrow(group_info_out),
    sum(group_info_out$group == "mild"), sum(group_info_out$group == "severe"),
    abundance_file, group_file))
}
