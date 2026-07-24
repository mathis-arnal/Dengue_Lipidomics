##########################################################################################
# THREE-GROUP (HEALTHY/MILD/SEVERE) BOXPLOTS FOR MILD-vs-SEVERE-SIGNIFICANT LIPIDS
##########################################################################################
#
# PBMC_04's boxplots only show mild vs severe (the two groups actually tested), since
# that script's processed_se never includes healthy samples. This script instead builds
# the three-group (healthy/mild/severe) processed_se -- the same one PBMC_01 uses -- and
# plots the same 11 lipids nominally significant in the Mild vs Severe comparison
# (analysis/PBMC/PostHoc_Mild_vs_Severe/D0/02.DiffExp/01.t-test_results/t-test_results.tsv),
# so healthy can be visually compared alongside mild and severe for each one. The
# significance bracket shown is still the Mild-vs-Severe t-test result specifically (the
# comparison that was actually tested) -- healthy is shown for visual context only, not
# as a third tested contrast.
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

base_dir <- "E:/Dengue_lipidomics/analysis/PBMC/PostHoc_Mild_vs_Severe/D0/02.DiffExp/03.Individual_lipid_boxplots_with_healthy"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("THREE-GROUP BOXPLOTS FOR MILD-vs-SEVERE-SIGNIFICANT LIPIDS\n")
cat("========================================\n")

loaded <- load_and_match_samples(
  "data/PBMC/group_information_table_healthy_vs_sick_patients_D0.tsv",
  "data/PBMC/PBMC_Lipid_abundance_data_D0.tsv"
)
abundance <- loaded$abundance
group_info <- loaded$group_info

annotated <- clean_and_annotate_lipids(abundance)

se_result <- build_processed_se(
  annotated$abundance_filtered, annotated$goslin_annotation, group_info,
  se_type = "de_multiple", exclude_missing_pct = 70, paired_sample = NULL
)
processed_se <- se_result$processed_se

abundance_data <- SummarizedExperiment::assay(processed_se)
col_data <- SummarizedExperiment::colData(processed_se)
cat("  Samples in three-group SE:", ncol(abundance_data), " (",
    paste(names(table(col_data$group)), table(col_data$group), sep = "=", collapse = ", "), ")\n")

ms_results <- read.delim("analysis/PBMC/PostHoc_Mild_vs_Severe/D0/02.DiffExp/01.t-test_results/t-test_results.tsv") %>%
  dplyr::arrange(pval)
cat("  Mild-vs-Severe nominally significant lipids:", nrow(ms_results), "\n")

group_colors <- c(healthy = "#4DAF4A", mild = "#377EB8", severe = "#E41A1C")
group_levels <- c("healthy", "mild", "severe")

for (i in seq_len(nrow(ms_results))) {
  lipid_name <- ms_results$feature[i]
  pval <- ms_results$pval[i]
  padj <- ms_results$padj[i]

  if (!(lipid_name %in% rownames(abundance_data))) {
    cat("  Skipping", lipid_name, "- not found in three-group processed data\n")
    next
  }

  pval_text <- if (pval < 0.001) "p < 0.001" else sprintf("p = %.4f", pval)
  sig_label <- if (pval < 0.001) "***" else if (pval < 0.01) "**" else if (pval < 0.05) "*" else "ns"

  plot_data <- data.frame(
    sample = colnames(abundance_data),
    abundance = as.numeric(abundance_data[lipid_name, ]),
    group = factor(col_data$group, levels = group_levels)
  )
  plot_data <- plot_data[!is.na(plot_data$abundance), ]

  p <- ggplot(plot_data, aes(x = group, y = abundance, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    geom_signif(
      comparisons = list(c("mild", "severe")),
      annotations = sig_label,
      textsize = 6, vjust = 0.5, map_signif_level = FALSE
    ) +
    labs(
      title = lipid_name,
      subtitle = paste0("Mild vs Severe: ", pval_text, ", padj = ", sprintf("%.3f", padj), " (healthy shown for context, not tested)"),
      x = "Group", y = "Log10 Abundance"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "none"
    ) +
    scale_fill_manual(values = group_colors)

  safe_name <- gsub("[^A-Za-z0-9_-]", "_", lipid_name)
  ggsave(
    filename = file.path(base_dir, paste0(sprintf("%03d", i), "_", safe_name, "_3group.png")),
    plot = p, width = 8, height = 6, dpi = 120
  )
  cat(sprintf("  [%d/%d] %s saved\n", i, nrow(ms_results), lipid_name))
}

cat("\nBoxplots saved to:", base_dir, "\n")
