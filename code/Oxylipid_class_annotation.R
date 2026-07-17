# Oxylipid class annotation
#
# Oxylipins (HETEs, HODEs, prostaglandins, resolvins, etc.) are fatty-acid oxidation
# products, not glycerophospholipids/sphingolipids, so rgoslin::parseLipidNames() can't
# assign them a real class/category. Feeding a fabricated annotation table into
# LipidSigR's as_summarized_experiment() was tried and does not work: it unconditionally
# runs lipid_annotation(), which performs real LION/ChEBI/LIPID MAPS ontology lookups and
# nomenclature parsing that crash on non-rgoslin names (confirmed empirically: row-count
# mismatch / 'x@assays' is not parallel to 'x'). So this table is used directly in plain
# R (Oxylipid_01_posthoc_mild_vs_severe.R, Oxylipid_02_longitudinal_analysis.R), not fed
# into any LipidSigR SummarizedExperiment function.
#
# `family` = parent PUFA the oxylipin is derived from, taken directly from the
# platform's own "Precursor" column in the Results sheet of the raw data
# (data/lipidsig_datasets/Oxylipines_raw/oxylipin_precursor_by_feature.tsv, written by
# Oxylipid_lipidomics_pretreatment.R) rather than inferred by hand — analogous to
# rgoslin's Lipid.Maps.Category.
# `class` = compound class (HETE, HODE, PG, LT, LX, RvD, RvE, PD, MaR, ...), hand-assigned
# from standard oxylipin/eicosanoid nomenclature since the platform doesn't provide this
# level of grouping directly — analogous to Lipid.Maps.Main.Class.

library(readr)
library(dplyr)
library(tibble)

out_dir <- "data/lipidsig_datasets/Oxylipines_raw"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

precursor_by_feature <- read_tsv(file.path(out_dir, "oxylipin_precursor_by_feature.tsv"), show_col_types = FALSE)

# Matched by stable ASCII prefix rather than the full string, since the omega (ꙍ6/ꙍ3)
# character in the source spreadsheet doesn't round-trip reliably through re-typed
# Unicode literals.
family_from_precursor <- function(precursor) {
  case_when(
    startsWith(precursor, "Linoleic acid")             ~ "LA-derived",
    startsWith(precursor, "Dihomo-gamma-linoleic acid") ~ "DGLA-derived",
    startsWith(precursor, "Arachidonic acid")           ~ "AA-derived",
    startsWith(precursor, "Alpha-linolenic acid")       ~ "ALA-derived",
    startsWith(precursor, "Eicosapentaenoic acid")      ~ "EPA-derived",
    startsWith(precursor, "Docosahexaenoic acid")       ~ "DHA-derived",
    TRUE ~ NA_character_
  )
}

class_lookup <- tribble(
  ~feature,          ~class,
  "9-HODE",          "HODE",
  "10-HODE",         "HODE",
  "13-HODE",         "HODE",
  "9,10-DiHOME",     "DiHOME",
  "12,13-DiHOME",    "DiHOME",
  "9-oxo-ODE",       "OxoODE",
  "13-oxo-ODE",      "OxoODE",
  "9,10,13-TriHOME", "TriHOME",
  "9,12,13-TriHOME", "TriHOME",
  "PGA1",            "PG",
  "5,6-EET",         "EET",
  "8,9-EET",         "EET",
  "11,12-EET",       "EET",
  "14,15-EET",       "EET",
  "5-HETE",          "HETE",
  "8-HETE",          "HETE",
  "12-HETE",         "HETE",
  "15-HETE",         "HETE",
  "5-oxo-ETE",       "OxoETE",
  "LTB4",            "LT",
  "LxA4",            "LX",
  "LXB4",            "LX",
  "8isoPGA2",        "IsoP",
  "PGD2",            "PG",
  "PGE2",            "PG",
  "PGF2a",           "PG",
  "11B-PGF2a",       "PG",
  "PGFM",            "PG",
  "15dPGJ2",         "PG",
  "6kPGF1a",         "PGI2",
  "TXB2",            "TX",
  "9-HOTrE",         "HOTrE",
  "13-HOTrE",        "HOTrE",
  "18-HEPE",         "HEPE",
  "LTB5",            "LT",
  "PGE3",            "PG",
  "RvE1",            "RvE",
  "5,6-DiHETE",      "DiHETE",
  "14-HDHA",         "HDHA",
  "17-HDHA",         "HDHA",
  "7(S)-MaR1",       "MaR",
  "PDx",             "PD",
  "PD1",             "PD",
  "RVD1",            "RvD",
  "RvD2",            "RvD",
  "RvD3",            "RvD",
  "RvD5",            "RvD"
)

oxylipin_class_annotation <- precursor_by_feature %>%
  mutate(family = family_from_precursor(precursor)) %>%
  select(feature, family) %>%
  left_join(class_lookup, by = "feature")

stopifnot(!anyNA(oxylipin_class_annotation$family), !anyNA(oxylipin_class_annotation$class))

write_tsv(oxylipin_class_annotation, file.path(out_dir, "oxylipin_class_annotation.tsv"))
cat(sprintf(
  "Saved %d-feature class annotation (%d families, %d classes) -> %s\n",
  nrow(oxylipin_class_annotation),
  n_distinct(oxylipin_class_annotation$family),
  n_distinct(oxylipin_class_annotation$class),
  file.path(out_dir, "oxylipin_class_annotation.tsv")))
