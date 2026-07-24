##########################################################################################
# PQN-NORMALIZED BOXPLOTS: ALL NOMINALLY SIGNIFICANT LIPIDS, ALL-PATIENTS HEALTHY vs
# DISEASED
##########################################################################################
#
# The PQN normalization check (see PBMC_06c_all_patients_check.R /
# analysis/PBMC/PostHoc_Healthy_vs_Diseased/D0/All_patients_check/) found FA 18:1 and
# PC 32:1 drop just below FDR significance under PQN (padj ~0.053) despite reaching it
# under Percentage normalization (padj = 0.042) -- both checks only ever produced a
# comparison TABLE of p-values, never actual boxplots of the PQN-normalized abundance
# values. This script builds the PQN-normalized SummarizedExperiment for the same
# 57-patient all-patients cohort and plots EVERY nominally significant lipid (pval<0.05,
# not just the two FDR-adjacent ones), sorted by significance -- matching the "nominal,
# not just FDR" boxplot convention already applied to every other comparison in this
# report -- using the same ggplot2 boxplot + jitter + significance-star style throughout.
##########################################################################################

suppressPackageStartupMessages({
  library(LipidSigR)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(rgoslin)
  library(ggplot2)
  library(ggsignif)
})

source("code/lipidomics_functions.R")

setwd("E:/Dengue_lipidomics")

base_dir <- "E:/Dengue_lipidomics/analysis/PBMC/PostHoc_Healthy_vs_Diseased_AllPatients/D0/02.DiffExp/03.Individual_lipid_boxplots_PQN"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("PQN BOXPLOTS: ALL NOMINALLY SIGNIFICANT LIPIDS, ALL-PATIENTS HEALTHY vs DISEASED\n")
cat("========================================\n")

loaded <- load_and_match_samples(
  "data/PBMC/group_information_table_healthy_vs_diseased_D0_all.tsv",
  "data/PBMC/PBMC_Lipid_abundance_data_D0.tsv"
)
abundance <- loaded$abundance
group_info <- loaded$group_info

annotated <- clean_and_annotate_lipids(abundance)

se_result <- build_processed_se(
  annotated$abundance_filtered, annotated$goslin_annotation, group_info,
  se_type = "de_two", exclude_missing_pct = 70, paired_sample = FALSE,
  normalization = "PQN"
)
se <- se_result$se
processed_se <- se_result$processed_se

# Data quality before/after PQN normalization -- same diagnostic pair shown for every
# other comparison in this report (Percentage-normalized there), just with PQN here.
cat("Assessing data quality (PQN)...\n")
quality_dir <- "E:/Dengue_lipidomics/analysis/PBMC/PostHoc_Healthy_vs_Diseased_AllPatients/D0/01.Profiling/00.Data_quality_PQN"
dir.create(quality_dir, recursive = TRUE, showWarnings = FALSE)
plot_quality <- plot_data_process(se, processed_se)

save_static_png(plot_quality$static_boxPlot_after,
  file.path(quality_dir, "Data_quality_boxplot_after_PQN.png"),
  width = 1200, height = 800
)
save_static_png(plot_quality$static_densityPlot_before,
  file.path(quality_dir, "Data_quality_density_before.png"),
  width = 1200, height = 800
)
save_static_png(plot_quality$static_densityPlot_after,
  file.path(quality_dir, "Data_quality_density_after_PQN.png"),
  width = 1200, height = 800
)
cat("  Data quality plots saved to:", quality_dir, "\n")

deSp_se <- deSp_twoGroup(
  processed_se,
  ref_group = "healthy", test = "t-test",
  significant = "pval", p_cutoff = 0.05, FC_cutoff = 1, transform = "log10"
)
res_list <- extract_summarized_experiment(deSp_se)
all_results <- as.data.frame(res_list$all_deSp_result)

abundance_data <- SummarizedExperiment::assay(processed_se)
col_data <- SummarizedExperiment::colData(processed_se)

sig_results <- all_results %>%
  dplyr::filter(pval < 0.05) %>%
  dplyr::arrange(pval)
target_lipids <- sig_results$feature
cat("  Nominally significant under PQN:", length(target_lipids), "->", paste(target_lipids, collapse = ", "), "\n")

for (i in seq_along(target_lipids)) {
  lipid_name <- target_lipids[i]
  pval <- all_results$pval[all_results$feature == lipid_name]
  padj <- all_results$padj[all_results$feature == lipid_name]
  cat(sprintf("  %s: pval=%.4g, padj=%.4g\n", lipid_name, pval, padj))

  pval_text <- if (pval < 0.001) "p < 0.001" else sprintf("p = %.4f (PQN)", pval)
  sig_label <- if (padj < 0.05) {
    "* (FDR-sig.)"
  } else if (pval < 0.001) {
    "***"
  } else if (pval < 0.01) {
    "**"
  } else if (pval < 0.05) {
    "*"
  } else {
    "ns"
  }

  plot_data <- data.frame(
    sample = colnames(abundance_data),
    abundance = as.numeric(abundance_data[lipid_name, ]),
    group = col_data$group
  )
  plot_data <- plot_data[!is.na(plot_data$abundance), ]

  p <- ggplot(plot_data, aes(x = group, y = abundance, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    geom_signif(
      comparisons = list(c("healthy", "diseased")),
      annotations = sig_label,
      textsize = 6, vjust = 0.5, map_signif_level = FALSE
    ) +
    labs(
      title = paste0(lipid_name, " (PQN-normalized)"),
      subtitle = paste0(pval_text, ", padj = ", sprintf("%.3f", padj)),
      x = "Group", y = "Log10 Abundance (PQN)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    scale_fill_manual(values = c("healthy" = "#4DAF4A", "diseased" = "#E41A1C"))

  safe_name <- gsub("[^A-Za-z0-9_-]", "_", lipid_name)
  ggsave(
    filename = file.path(base_dir, paste0(sprintf("%03d", i), "_", safe_name, "_PQN.png")),
    plot = p, width = 8, height = 6, dpi = 120
  )
}

cat("\nBoxplots saved to:", base_dir, "\n")
