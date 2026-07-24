##########################################################################################
# PLASMA POST-HOC: D0 vs D60 (CONVALESCENCE AS CONTROL), WITHIN MILD AND WITHIN SEVERE
##########################################################################################
#
# Requested addition to the plasma analysis: since there is no longitudinal healthy-
# control arm in this project (see 05_longitudinal_analysis.R header), D60 (~2 months
# post-enrollment, well past acute illness) is used here as a within-patient "recovered"
# baseline instead. This directly asks "how far has each patient's lipidome shifted from
# its own later, presumably-recovered state, at the time of acute presentation" -- run
# SEPARATELY for mild and severe patients (not pooled), since that's what was asked and
# since a pooled test would conflate two different-severity patient populations' distinct
# baselines.
#
# Statistically this is a PAIRED two-group comparison (same patient, D0 vs D60), unlike
# the cross-sectional Healthy/Mild/Severe posthoc scripts (02/03/04), which compare
# different, unrelated patients. LipidSigR's deSp_twoGroup(paired_sample=TRUE) is used
# accordingly, with patient_id as the pairing key -- this is a stronger, better-powered
# design than an unpaired test for the same n, since each patient serves as their own
# control.
#
# DATA SOURCE: same as 05_longitudinal_analysis.R --
# data/lipidsig_datasets/Longitudinal_matched/Lipid_abundance_data_{D0,D60}.tsv +
# group_information_table_{D0,D60}.xlsx (NOT the old Lipid_abundance_data_D0_vs_D60_severe.tsv
# / group_information_table_severe_D0_D60.xlsx files that used to live alongside these in
# data/lipidsig_datasets/Longitudinal_matched/ -- their group_info sample_name/label_name/group
# columns didn't line up with each other, likely an abandoned draft; deleted 2026-07-24).
##########################################################################################

suppressPackageStartupMessages({
  library(LipidSigR)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(readxl)
  library(rgoslin)
  library(stringr)
  library(ggplot2)
  library(ggsignif)
})

source("code/lipidomics_functions.R")

setwd("E:/Dengue_lipidomics")

cat("\n========================================\n")
cat("PLASMA POST-HOC: D0 vs D60, BY SEVERITY\n")
cat("========================================\n")


##########################################################################################
# 1. DATA LOADING (D0 + D60 ONLY) -- mirrors 05_longitudinal_analysis.R Section 1
##########################################################################################

cat("\n=== 1. DATA LOADING ===\n")

d0 <- read_tsv("data/lipidsig_datasets/Longitudinal_matched/Lipid_abundance_data_D0.tsv", show_col_types = FALSE)
d60 <- read_tsv("data/lipidsig_datasets/Longitudinal_matched/Lipid_abundance_data_D60.tsv", show_col_types = FALSE)
stopifnot(identical(d0$feature, d60$feature))

abund_combined <- dplyr::full_join(d0, d60, by = "feature")
cat("  Combined abundance table:", nrow(abund_combined), "features x", ncol(abund_combined) - 1, "sample-timepoints\n")

group_combined <- bind_rows(
  read_excel("data/lipidsig_datasets/Longitudinal_matched/group_information_table_D0.xlsx"),
  read_excel("data/lipidsig_datasets/Longitudinal_matched/group_information_table_D60.xlsx")
)

# Same constant-value cleaning as 05_longitudinal_analysis.R: 2 patients (JV-048, KT-565)
# have an all-zero/constant D60 column (no real D60 draw), which as_summarized_experiment()
# cannot handle cleanly when combined with any constant feature.
mat <- as.matrix(abund_combined[, -1])
rownames(mat) <- abund_combined$feature
feat_const <- apply(mat, 1, function(x) length(unique(x)) <= 1)
samp_const <- apply(mat, 2, function(x) length(unique(x)) <= 1)
if (any(feat_const)) cat("  Dropping", sum(feat_const), "constant-value feature(s)\n")
if (any(samp_const)) cat("  Dropping", sum(samp_const), "constant-value sample-timepoint(s):", paste(colnames(mat)[samp_const], collapse = ", "), "\n")

abund_clean <- abund_combined %>% dplyr::filter(!feat_const)
keep_samples <- colnames(mat)[!samp_const]
abund_clean <- abund_clean[, c("feature", keep_samples)]
group_clean <- group_combined %>% dplyr::filter(sample_name %in% keep_samples)

sample_meta <- group_clean %>%
  mutate(
    timepoint = case_when(
      str_detect(sample_name, "D60$") ~ "D60",
      str_detect(sample_name, "D0$")  ~ "D0",
      TRUE ~ NA_character_
    ),
    patient_id = str_remove(sample_name, "D(60|0)$")
  ) %>%
  rename(severity = group) %>%
  select(sample_name, patient_id, timepoint, severity)

# Keep only patients with BOTH D0 and D60 (required for a paired test)
paired_patients <- sample_meta %>%
  count(patient_id) %>%
  filter(n == 2) %>%
  pull(patient_id)
cat("  Patients with both D0 and D60:", length(paired_patients), "\n")
sample_meta <- sample_meta %>% filter(patient_id %in% paired_patients)
cat("  By severity:\n")
print(sample_meta %>% distinct(patient_id, severity) %>% count(severity))


##########################################################################################
# 2. LIPID ANNOTATION (rgoslin) -- done once, shared across both severity runs
##########################################################################################

cat("\n=== 2. LIPID ANNOTATION ===\n")

parse_lipid <- rgoslin::parseLipidNames(lipidNames = abund_clean$feature)
recognized_lipids <- parse_lipid$Original.Name[which(parse_lipid$Grammar != "NOT_PARSEABLE")]
abundance_filtered_all <- abund_clean %>% dplyr::filter(feature %in% recognized_lipids)
goslin_annotation <- parse_lipid %>% dplyr::filter(Original.Name %in% recognized_lipids)
cat("  Recognized lipids:", length(recognized_lipids), "/", nrow(abund_clean), "\n")


##########################################################################################
# 3. PER-SEVERITY PAIRED D0 vs D60 ANALYSIS
##########################################################################################

run_severity_analysis <- function(severity_level) {
  cat("\n\n########################################\n")
  cat("SEVERITY:", toupper(severity_level), "\n")
  cat("########################################\n")

  base_dir <- file.path("E:/Dengue_lipidomics/analysis", paste0("D0_vs_D60_", str_to_title(severity_level)))
  create_output_dirs(base_dir, c(
    "01.Profiling/00.Data_quality",
    "01.Profiling/02.Dimensionality_reduction/PCA",
    "01.Profiling/03.Correlation_Heatmap",
    "01.Profiling/04.Lipid_characteristics",
    "02.DiffExp/01.t-test_results",
    "02.DiffExp/02.Visualizations/Volcano",
    "02.DiffExp/02.Visualizations/MA_plot",
    "02.DiffExp/03.Individual_lipid_boxplots",
    "03.Enrichment/LSEA",
    "03.Enrichment/ORA"
  ))

  meta_sev <- sample_meta %>% filter(severity == severity_level)
  keep_cols <- meta_sev$sample_name
  cat("  Samples (", length(keep_cols), "= ", n_distinct(meta_sev$patient_id), "patients x 2 timepoints):\n")

  abundance_sev <- abundance_filtered_all[, c("feature", keep_cols)]

  # as_summarized_experiment()'s paired_sample=TRUE requires the 'pair' column to be
  # sequential integers 1..N (one per patient), not patient ID strings.
  group_info <- meta_sev %>%
    transmute(sample_name, label_name = sample_name, group = timepoint,
              pair = as.integer(factor(patient_id, levels = unique(patient_id)))) %>%
    arrange(match(sample_name, keep_cols))
  # as_summarized_experiment() needs group_info rows in the same order as abundance columns
  group_info <- group_info[match(colnames(abundance_sev)[-1], group_info$sample_name), ]
  stopifnot(identical(group_info$sample_name, colnames(abundance_sev)[-1]))

  se_result <- build_processed_se(
    abundance_sev, goslin_annotation, group_info,
    se_type = "de_two", exclude_missing_pct = 70, paired_sample = TRUE
  )
  se <- se_result$se
  processed_se <- se_result$processed_se

  # --- Profiling ---
  cat("  Assessing data quality...\n")
  plot_quality <- plot_data_process(se, processed_se)
  save_static_png(plot_quality$static_boxPlot_before, file.path(base_dir, "01.Profiling/00.Data_quality/Data_quality_boxplot_before.png"), width = 1200, height = 800)
  save_static_png(plot_quality$static_boxPlot_after, file.path(base_dir, "01.Profiling/00.Data_quality/Data_quality_boxplot_after.png"), width = 1200, height = 800)
  save_static_png(plot_quality$static_densityPlot_before, file.path(base_dir, "01.Profiling/00.Data_quality/Data_quality_density_before.png"), width = 1200, height = 800)
  save_static_png(plot_quality$static_densityPlot_after, file.path(base_dir, "01.Profiling/00.Data_quality/Data_quality_density_after.png"), width = 1200, height = 800)

  cat("  Performing PCA...\n")
  plot_pca <- dr_pca(processed_se, scaling = TRUE, centering = TRUE, clustering = "kmeans", cluster_num = 2)
  pca_timepoint <- plot_pca_by_severity(plot_pca, group_info)
  save_static_png(pca_timepoint$plot, file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/PCA_plot_by_timepoint.png"), width = 1000, height = 800)
  save_interactive_rds(pca_timepoint$plot_interactive, file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/PCA_plot_by_timepoint_interactive.rds"))
  save_static_png(plot_pca$static_pca, file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/PCA_plot.png"), width = 1000, height = 800)
  save_interactive_rds(plot_pca$interactive_pca, file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/PCA_plot_interactive.rds"))
  write_tsv(plot_pca$pca_rotated_data, file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/table_pca_rotated_data.tsv"))

  cat("  Computing correlation heatmap...\n")
  plot_corr <- heatmap_correlation(processed_se, char = NULL, transform = "log10", correlation = "spearman", distfun = "euclidean", hclustfun = "ward.D2", type = "sample")
  save_static_png(plot_corr$static_heatmap, file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/Correlation_heatmap_samples.png"), width = 1200, height = 1000)
  save_interactive_rds(plot_corr$interactive_heatmap, file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/Correlation_heatmap_samples_interactive.rds"))

  cat("  Lipid characteristics...\n")
  result_class <- lipid_profiling(processed_se, char = "class")
  save_static_png(result_class$static_lipid_composition, file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Lipid_composition_per_sample.png"), width = 1000, height = 1000)

  # --- Differential expression (paired t-test, D60 = reference/"control") ---
  cat("  Running paired t-test: D0 vs D60 (D60 = reference)...\n")
  deSp_se <- deSp_twoGroup(
    processed_se,
    ref_group = "D60", test = "t-test",
    significant = "pval", p_cutoff = 0.05, FC_cutoff = 1, transform = "log10"
  )
  res_list <- extract_summarized_experiment(deSp_se)
  desp_df <- as.data.frame(res_list$sig_deSp_result)
  cat("  Nominally significant (p<0.05):", nrow(desp_df), "\n")
  cat("  FDR-significant (padj<0.05):", sum(desp_df$padj < 0.05, na.rm = TRUE), "\n")
  write_tsv(desp_df, file.path(base_dir, "02.DiffExp/01.t-test_results/t-test_results.tsv"))

  all_df <- as.data.frame(res_list$all_deSp_result)
  write_tsv(all_df, file.path(base_dir, "02.DiffExp/01.t-test_results/t-test_all_results.tsv"))

  deSp_plot <- plot_deSp_twoGroup(deSp_se)
  save_static_png(deSp_plot$static_volcanoPlot, file.path(base_dir, "02.DiffExp/02.Visualizations/Volcano/Volcano_plot.png"))
  save_interactive_rds(deSp_plot$interactive_volcanoPlot, file.path(base_dir, "02.DiffExp/02.Visualizations/Volcano/Volcano_plot_interactive.rds"))
  save_static_png(deSp_plot$static_maPlot, file.path(base_dir, "02.DiffExp/02.Visualizations/MA_plot/MA_plot.png"))
  save_interactive_rds(deSp_plot$interactive_maPlot, file.path(base_dir, "02.DiffExp/02.Visualizations/MA_plot/MA_plot_interactive.rds"))

  # --- Boxplots: all nominally significant lipids, sorted by significance, paired lines
  # connecting each patient's D0->D60 values (more informative than an unpaired boxplot
  # for a paired design) ---
  if (nrow(desp_df) > 0) {
    sig_sorted <- desp_df %>% arrange(pval)
    abundance_data <- SummarizedExperiment::assay(processed_se)
    col_data <- SummarizedExperiment::colData(processed_se)
    cat("  Creating paired boxplots for", nrow(sig_sorted), "significant lipids...\n")

    for (i in seq_len(nrow(sig_sorted))) {
      lipid_name <- sig_sorted$feature[i]
      pval <- sig_sorted$pval[i]
      padj <- sig_sorted$padj[i]
      pval_text <- if (pval < 0.001) "p < 0.001" else sprintf("p = %.4f", pval)
      sig_label <- if (padj < 0.05) "* (FDR-sig.)" else if (pval < 0.001) "***" else if (pval < 0.01) "**" else "*"

      plot_data <- data.frame(
        sample = colnames(abundance_data),
        abundance = as.numeric(abundance_data[lipid_name, ]),
        timepoint = factor(col_data$group, levels = c("D0", "D60")),
        patient_id = col_data$pair
      )
      plot_data <- plot_data[!is.na(plot_data$abundance), ]

      p <- ggplot(plot_data, aes(x = timepoint, y = abundance)) +
        geom_line(aes(group = patient_id), alpha = 0.3, color = "grey50") +
        geom_boxplot(aes(fill = timepoint), outlier.shape = NA, alpha = 0.7, width = 0.5) +
        geom_jitter(aes(color = timepoint), width = 0.05, alpha = 0.7, size = 2) +
        geom_signif(
          comparisons = list(c("D0", "D60")), annotations = sig_label,
          textsize = 6, vjust = 0.5, map_signif_level = FALSE
        ) +
        labs(
          title = lipid_name,
          subtitle = paste0("D0 vs D60 (paired), ", pval_text, ", padj = ", sprintf("%.3f", padj)),
          x = "Timepoint", y = "Log10 Abundance"
        ) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 9), legend.position = "none") +
        scale_fill_manual(values = c(D0 = "#E41A1C", D60 = "#4DAF4A")) +
        scale_color_manual(values = c(D0 = "#E41A1C", D60 = "#4DAF4A"))

      safe_name <- gsub("[^A-Za-z0-9_-]", "_", lipid_name)
      ggsave(
        filename = file.path(base_dir, "02.DiffExp/03.Individual_lipid_boxplots", paste0(sprintf("%03d", i), "_", safe_name, ".png")),
        plot = p, width = 8, height = 6, dpi = 120
      )
    }
  }

  # --- Enrichment ---
  run_lsea_ora(
    deSp_se,
    char_list = list(by_class = "class", by_category = "Category"),
    lsea_dir = file.path(base_dir, "03.Enrichment/LSEA"),
    ora_dir = file.path(base_dir, "03.Enrichment/ORA")
  )

  cat("\n  ✓", severity_level, "analysis complete. Output:", base_dir, "\n")
  invisible(list(severity = severity_level, n_sig = nrow(desp_df), n_fdr = sum(desp_df$padj < 0.05, na.rm = TRUE)))
}

results <- lapply(c("mild", "severe"), run_severity_analysis)

cat("\n\n========================================\n")
cat("D0 vs D60 ANALYSIS COMPLETE (BOTH SEVERITIES)\n")
cat("========================================\n")
for (r in results) cat(" ", r$severity, ": ", r$n_sig, "nominally significant,", r$n_fdr, "FDR-significant\n")
