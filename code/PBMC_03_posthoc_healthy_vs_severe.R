##########################################################################################
# LIPIDOMICS ANALYSIS: POST-HOC HEALTHY vs SEVERE
##########################################################################################
#
# COMPLETE ANALYSIS WORKFLOW:
# 1. Data Loading & Preparation
# 2. Profiling (Data Quality, Cross-Sample Variability, Dimensionality Reduction)
# 3. Correlation Analysis
# 4. Lipid Characteristics
# 5. Differential Expression Analysis (t-test)
# 6. Individual Lipid Boxplots
# 7. Lipid Set Enrichment Analysis (LSEA)
# 8. Over Representation Analysis (ORA)
#
# This script performs pairwise comparison: Healthy vs Severe
##########################################################################################

# Load required libraries
library(LipidSigR)
library(SummarizedExperiment)
library(dplyr)
library(readr)
library(rgoslin)
library(ggplot2)
library(ggsignif)
library(mixOmics)

source("code/lipidomics_functions.R")

# Set working directory
setwd("E:/Dengue_lipidomics")

# Define base directory for outputs
base_dir <- "E:/Dengue_lipidomics/analysis/PBMC/PostHoc_Healthy_vs_Severe/D0"

cat("\n========================================\n")
cat("POST-HOC ANALYSIS: HEALTHY vs SEVERE\n")
cat("========================================\n")


##########################################################################################
# 1. DATA LOADING & PREPARATION
##########################################################################################

cat("\n=== 1. DATA LOADING & PREPARATION ===\n")

create_output_dirs(base_dir, c(
  "01.Profiling/00.Data_quality",
  "01.Profiling/01.Cross-sample_variability",
  "01.Profiling/02.Dimensionality_reduction/PCA",
  "01.Profiling/02.Dimensionality_reduction/t-SNE",
  "01.Profiling/02.Dimensionality_reduction/UMAP",
  "01.Profiling/03.Correlation_Heatmap/by_samples",
  "01.Profiling/03.Correlation_Heatmap/by_category",
  "01.Profiling/03.Correlation_Heatmap/by_class",
  "01.Profiling/04.Lipid_characteristics",
  "02.DiffExp/01.t-test_results",
  "02.DiffExp/02.Visualizations/Volcano",
  "02.DiffExp/02.Visualizations/MA_plot",
  "02.DiffExp/02.Visualizations/Heatmap",
  "02.DiffExp/02.Visualizations/PLS-DA",
  "02.DiffExp/03.Individual_lipid_boxplots",
  "03.Enrichment/LSEA",
  "03.Enrichment/ORA"
))

cat("Loading group information and abundance data...\n")
loaded <- load_and_match_samples(
  "data/PBMC/group_information_table_healthy_vs_SEVERE_D0.tsv",
  "data/PBMC/PBMC_Lipid_abundance_data_D0.tsv"
)
abundance <- loaded$abundance
group_info <- loaded$group_info

annotated <- clean_and_annotate_lipids(abundance)

# Construct SummarizedExperiment object + process
se_result <- build_processed_se(
  annotated$abundance_filtered, annotated$goslin_annotation, group_info,
  se_type = "de_two", exclude_missing_pct = 50, paired_sample = FALSE
)
se <- se_result$se
processed_se <- se_result$processed_se


##########################################################################################
# 2. PROFILING
##########################################################################################

cat("\n=== 2. PROFILING ===\n")

# Data quality assessment - compare before/after processing
cat("Assessing data quality...\n")
plot_quality <- plot_data_process(se, processed_se)

save_static_png(plot_quality$static_boxPlot_before,
  file.path(base_dir, "01.Profiling/00.Data_quality/Data_quality_boxplot_before.png"),
  width = 1200, height = 800
)
save_static_png(plot_quality$static_boxPlot_after,
  file.path(base_dir, "01.Profiling/00.Data_quality/Data_quality_boxplot_after.png"),
  width = 1200, height = 800
)
save_static_png(plot_quality$static_densityPlot_before,
  file.path(base_dir, "01.Profiling/00.Data_quality/Data_quality_density_before.png"),
  width = 1200, height = 800
)
save_static_png(plot_quality$static_densityPlot_after,
  file.path(base_dir, "01.Profiling/00.Data_quality/Data_quality_density_after.png"),
  width = 1200, height = 800
)

# Cross-sample variability
cat("Analyzing cross-sample variability...\n")
plot_variability <- cross_sample_variability(processed_se)

save_static_png(plot_variability$static_lipid_number_barPlot,
  file.path(base_dir, "01.Profiling/01.Cross-sample_variability/Expressed_lipid_numbers.png"),
  width = 1200, height = 800
)
save_static_png(plot_variability$static_lipid_distribution,
  file.path(base_dir, "01.Profiling/01.Cross-sample_variability/Lipid_abundance_distribution.png"),
  width = 1200, height = 800
)

# PCA
cat("Performing PCA...\n")
plot_pca <- dr_pca(
  processed_se,
  scaling = TRUE,
  centering = TRUE,
  clustering = "kmeans",
  cluster_num = 2
)

# Severity-colored PCA (color = true clinical group, shape = unsupervised k-means
# cluster, 95% ellipse per group) -- see plot_pca_by_severity() in
# lipidomics_functions.R for why clustering="group_info" can't be used directly here.
pca_severity <- plot_pca_by_severity(plot_pca, group_info)

save_static_png(pca_severity$plot,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/PCA_plot_clinical_group.png"),
  width = 1000, height = 800
)
save_interactive_rds(pca_severity$plot_interactive,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/PCA_plot_clinical_group_interactive.rds"))

save_static_png(plot_pca$static_pca,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/PCA_plot.png"),
  width = 1000, height = 800
)
save_static_png(plot_pca$static_screePlot,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Scree_plot.png"),
  width = 1000, height = 800
)

save_interactive_rds(plot_pca$interactive_pca,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/PCA_plot_interactive.rds"))
save_interactive_rds(plot_pca$interactive_screePlot,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Scree_plot_interactive.rds"))
save_interactive_rds(plot_pca$interactive_variablePlot,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Variable_correlation_interactive.rds"))

write_tsv(plot_pca$pca_rotated_data,
  file = file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/table_pca_rotated_data.tsv")
)

# t-SNE
cat("Performing t-SNE...\n")
plot_tsne <- dr_tsne(
  processed_se,
  pca = TRUE,
  perplexity = 5,
  clustering = "kmeans",
  cluster_num = 2
)

save_static_png(plot_tsne$static_tsne,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/t-SNE/tSNE_plot.png"),
  width = 1000, height = 800
)
save_interactive_rds(plot_tsne$interactive_tsne,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/t-SNE/t-SNE_plot_interactive.rds"))

write_tsv(plot_tsne$tsne_result,
  file = file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/t-SNE/table_tsne_results.tsv")
)

# UMAP
cat("Performing UMAP...\n")
plot_umap <- dr_umap(
  processed_se,
  n_neighbors = 15,
  scaling = TRUE,
  umap_metric = "euclidean",
  clustering = "kmeans",
  cluster_num = 2
)

save_static_png(plot_umap$static_umap,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/UMAP/UMAP_plot.png"),
  width = 1000, height = 800
)
save_interactive_rds(plot_umap$interactive_umap,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/UMAP/UMAP_plot_interactive.rds"))

write_tsv(plot_umap$umap_result,
  file = file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/UMAP/UMAP_results.tsv")
)

cat("  ✓ Profiling complete\n")


##########################################################################################
# 3. CORRELATION ANALYSIS
##########################################################################################

cat("\n=== 3. CORRELATION ANALYSIS ===\n")

# Correlation by samples
cat("Computing correlation by samples...\n")
plot_corr_samples <- heatmap_correlation(
  processed_se,
  char = NULL,
  transform = "log10",
  correlation = "spearman",
  distfun = "euclidean",
  hclustfun = "ward.D2",
  type = "sample"
)

save_static_png(plot_corr_samples$static_heatmap,
  file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_samples/Correlation_heatmap_samples.png"),
  width = 1200, height = 1000
)
save_interactive_rds(plot_corr_samples$interactive_heatmap,
  file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_samples/Correlation_heatmap_samples_interactive.rds"))

# Correlation by category
cat("Computing correlation by category...\n")
plot_corr_category <- heatmap_correlation(
  processed_se,
  char = "Category",
  transform = "log10",
  correlation = "spearman",
  distfun = "euclidean",
  hclustfun = "ward.D2",
  type = "class"
)

save_static_png(plot_corr_category$static_heatmap,
  file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_category/Correlation_heatmap_category.png"),
  width = 1200, height = 1000
)
save_interactive_rds(plot_corr_category$interactive_heatmap,
  file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_category/Correlation_heatmap_category_interactive.rds"))

# Correlation by class
cat("Computing correlation by class...\n")
plot_corr_class <- heatmap_correlation(
  processed_se,
  char = "class",
  transform = "log10",
  correlation = "spearman",
  distfun = "euclidean",
  hclustfun = "ward.D2",
  type = "class"
)

save_static_png(plot_corr_class$static_heatmap,
  file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_class/Correlation_heatmap_class.png"),
  width = 1500, height = 1200
)
save_interactive_rds(plot_corr_class$interactive_heatmap,
  file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_class/Correlation_heatmap_class_interactive.rds"))

cat("  ✓ Correlation analysis complete\n")


##########################################################################################
# 4. LIPID CHARACTERISTICS
##########################################################################################

cat("\n=== 4. LIPID CHARACTERISTICS ===\n")

# Characteristics by class
cat("Analyzing lipid characteristics by class...\n")
result_class <- lipid_profiling(processed_se, char = "class")

save_static_png(result_class$static_char_barPlot,
  file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Lipid_characteristics_by_class.png"))
save_static_png(result_class$static_lipid_composition,
  file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Lipid_composition_per_sample.png"),
  width = 1000, height = 1000
)

# Characteristics by Total.C
cat("Analyzing lipid characteristics by Total.C...\n")
result_totalC <- lipid_profiling(processed_se, char = "Total.C")

save_static_png(result_totalC$static_char_barPlot,
  file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Lipid_characteristics_by_TotalC.png"))

# Characteristics by Total.DB
cat("Analyzing lipid characteristics by Total.DB...\n")
result_totalDB <- lipid_profiling(processed_se, char = "Total.DB")

save_static_png(result_totalDB$static_char_barPlot,
  file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Lipid_characteristics_by_TotalDB.png"))

# Characteristics by Category
cat("Analyzing lipid characteristics by Category...\n")
result_category <- lipid_profiling(processed_se, char = "Category")

save_static_png(result_category$static_char_barPlot,
  file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Lipid_characteristics_by_Category.png"))

cat("  ✓ Lipid characteristics analysis complete\n")


##########################################################################################
# 5. DIFFERENTIAL EXPRESSION ANALYSIS (T-TEST)
##########################################################################################

cat("\n=== 5. DIFFERENTIAL EXPRESSION ANALYSIS (t-test) ===\n")

# Two-group differential expression using t-test
cat("Running t-test: Healthy vs Severe...\n")
deSp_se <- deSp_twoGroup(
  processed_se,
  ref_group = "healthy",
  test = "t-test",
  significant = "pval",
  p_cutoff = 0.05,
  FC_cutoff = 1,
  transform = "log10"
)

# Extract results
res_list <- extract_summarized_experiment(deSp_se)

# For two-group analysis, check multiple possible result structures
if (!is.null(res_list$sig_deSp_result)) {
  desp_df <- as.data.frame(res_list$sig_deSp_result)
} else if (!is.null(res_list$all_deSp_result)) {
  desp_df <- as.data.frame(res_list$all_deSp_result)
} else {
  cat("  Warning: Could not find expected result structure.\n")
  cat("  Available result fields:\n")
  print(names(res_list))
  desp_df <- NULL
}

if (!is.null(desp_df)) {
  cat("  Total lipids tested:", nrow(desp_df), "\n")
  cat("  Significant lipids (padj < 0.05):", sum(desp_df$padj < 0.05, na.rm = TRUE), "\n")
  cat("  Up-regulated in Severe:", sum(desp_df$padj < 0.05 & desp_df$log2FC > 0, na.rm = TRUE), "\n")
  cat("  Down-regulated in Severe:", sum(desp_df$padj < 0.05 & desp_df$log2FC < 0, na.rm = TRUE), "\n")

  # Save results
  write_tsv(desp_df,
    file = file.path(base_dir, "02.DiffExp/01.t-test_results/t-test_results.tsv")
  )
}

# Plot differential expression results
deSp_plot <- plot_deSp_twoGroup(deSp_se)

# Volcano plot
save_static_png(deSp_plot$static_volcanoPlot,
  file.path(base_dir, "02.DiffExp/02.Visualizations/Volcano/Volcano_plot.png"))
save_interactive_html_and_rds(
  deSp_plot$interactive_volcanoPlot,
  file.path(base_dir, "02.DiffExp/02.Visualizations/Volcano/Volcano_plot_interactive.html"),
  file.path(base_dir, "02.DiffExp/02.Visualizations/Volcano/Volcano_plot_interactive.rds")
)
if (!is.null(deSp_plot$interactive_volcanoPlot)) cat("  Interactive volcano plot saved\n")

# MA plot
save_static_png(deSp_plot$static_maPlot,
  file.path(base_dir, "02.DiffExp/02.Visualizations/MA_plot/MA_plot.png"))
save_interactive_html_and_rds(
  deSp_plot$interactive_maPlot,
  file.path(base_dir, "02.DiffExp/02.Visualizations/MA_plot/MA_plot_interactive.html"),
  file.path(base_dir, "02.DiffExp/02.Visualizations/MA_plot/MA_plot_interactive.rds")
)
if (!is.null(deSp_plot$interactive_maPlot)) cat("  Interactive MA plot saved\n")

# Lollipop chart
save_static_png(deSp_plot$static_de_lipid,
  file.path(base_dir, "02.DiffExp/02.Visualizations/Lollipop_chart.png"))

# Heatmap
save_static_png(deSp_plot$static_heatmap,
  file.path(base_dir, "02.DiffExp/02.Visualizations/Heatmap/Heatmap_significant_lipids.png"),
  width = 1500, height = 1200
)

cat("  ✓ Differential expression analysis complete\n")


##########################################################################################
# 5A. PLS-DA ANALYSIS
##########################################################################################

cat("\n=== 5A. PLS-DA ANALYSIS ===\n")

# Perform PLS-DA. dr_plsda() requires at least 3 padj-significant lipids (via
# LipidSigR's internal .deSig_data(sigN=3)) -- comparisons with very few significant
# hits (as Healthy vs Severe/Mild vs Severe can be) will fail here. Wrapped in
# tryCatch() so that data-dependent case degrades gracefully instead of halting the
# whole script.
cat("Running PLS-DA...\n")
result_plsda <- tryCatch(
  dr_plsda(
    deSp_se,
    ncomp = 2,
    scaling = TRUE,
    clustering = "group_info",
    cluster_num = 2,
    kmedoids_metric = NULL,
    distfun = NULL,
    hclustfun = NULL,
    eps = NULL,
    minPts = NULL
  ),
  error = function(e) {
    cat("  Warning: PLS-DA could not be run -", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(result_plsda)) {
  # Save PLS-DA results tables
  write.csv(result_plsda$plsda_result,
    file = file.path(base_dir, "02.DiffExp/02.Visualizations/PLS-DA/Table_of_sample_variate.csv"),
    row.names = FALSE
  )
  write.csv(result_plsda$table_plsda_loading,
    file = file.path(base_dir, "02.DiffExp/02.Visualizations/PLS-DA/Table_of_sample_loading.csv"),
    row.names = FALSE
  )

  save_static_png(result_plsda$static_plsda,
    file.path(base_dir, "02.DiffExp/02.Visualizations/PLS-DA/PLSDA_plot.png"),
    width = 1200, height = 800
  )
  save_interactive_rds(result_plsda$interacitve_plsda,
    file.path(base_dir, "02.DiffExp/02.Visualizations/PLS-DA/PLSDA_plot_interactive.rds"))

  save_static_png(result_plsda$static_loadingPlot,
    file.path(base_dir, "02.DiffExp/02.Visualizations/PLS-DA/Loading_plot.png"),
    width = 800, height = 800
  )
  save_interactive_rds(result_plsda$interactive_loadingPlot,
    file.path(base_dir, "02.DiffExp/02.Visualizations/PLS-DA/Loading_plot_interactive.rds"))

  cat("  ✓ PLS-DA analysis complete\n")
}


##########################################################################################
# 5B. DIFFERENTIAL EXPRESSION OF LIPID CHARACTERISTICS
##########################################################################################

cat("\n=== 5B. DIFFERENTIAL EXPRESSION OF LIPID CHARACTERISTICS ===\n")

cat("Analyzing differential expression of lipid characteristics...\n")

# By Total Carbon (Total.C)
cat("  - By Total Carbon\n")
deChar_totalC_se <- deChar_twoGroup(
  processed_se,
  char = "Total.C",
  ref_group = "healthy",
  test = "t-test",
  significant = "pval",
  p_cutoff = 0.05,
  FC_cutoff = 1,
  transform = "log10"
)
deChar_totalC_plot <- plot_deChar_twoGroup(deChar_totalC_se)
save_static_png(deChar_totalC_plot$static_char_barPlot,
  file.path(base_dir, "02.DiffExp/02.Visualizations/deChar_TotalC_barplot.png"))

# By Category
cat("  - By Category\n")
deChar_category_se <- deChar_twoGroup(
  processed_se,
  char = "Category",
  ref_group = "healthy",
  test = "t-test",
  significant = "pval",
  p_cutoff = 0.05,
  FC_cutoff = 1,
  transform = "log10"
)
deChar_category_plot <- plot_deChar_twoGroup(deChar_category_se)
save_static_png(deChar_category_plot$static_char_barPlot,
  file.path(base_dir, "02.DiffExp/02.Visualizations/deChar_Category_barplot.png"))

# By Class
cat("  - By Class\n")
deChar_class_se <- deChar_twoGroup(
  processed_se,
  char = "class",
  ref_group = "healthy",
  test = "t-test",
  significant = "pval",
  p_cutoff = 0.05,
  FC_cutoff = 1,
  transform = "log10"
)
deChar_class_plot <- plot_deChar_twoGroup(deChar_class_se)
save_static_png(deChar_class_plot$static_char_barPlot,
  file.path(base_dir, "02.DiffExp/02.Visualizations/deChar_Class_barplot.png"))

cat("  ✓ Differential expression of lipid characteristics complete\n")


##########################################################################################
# 6. INDIVIDUAL LIPID BOXPLOTS
##########################################################################################

cat("\n=== 6. INDIVIDUAL LIPID BOXPLOTS ===\n")

# Get significant lipids
if (exists("desp_df") && !is.null(desp_df) && nrow(desp_df) > 0) {
  # Check which p-value column is available
  # Nominal (uncorrected) p-value first -- boxplots should cover every nominally
  # significant lipid, not just the (usually much smaller/empty) FDR-significant set.
  pval_col <- if ("pval" %in% colnames(desp_df)) {
    "pval"
  } else if ("padj" %in% colnames(desp_df)) {
    "padj"
  } else {
    NULL
  }

  if (!is.null(pval_col)) {
    sig_lipids <- desp_df %>%
      filter(.data[[pval_col]] < 0.05) %>%
      arrange(.data[[pval_col]])

    n_sig <- nrow(sig_lipids)
    cat("  Creating boxplots for", n_sig, "significant lipids...\n")
    cat("  Using p-value column:", pval_col, "\n")

    if (n_sig > 0) {
      max_plots <- min(50, n_sig)
      cat("  (Showing top", max_plots, "most significant lipids)\n")

      # Get abundance data and group info from processed_se
      abundance_data <- SummarizedExperiment::assay(processed_se)
      col_data <- SummarizedExperiment::colData(processed_se)

      success_count <- 0
      failed_count <- 0

      for (i in 1:max_plots) {
        lipid_name <- if ("feature" %in% colnames(sig_lipids)) {
          sig_lipids$feature[i]
        } else if ("lipid" %in% colnames(sig_lipids)) {
          sig_lipids$lipid[i]
        } else {
          rownames(sig_lipids)[i]
        }

        safe_name <- gsub("[^A-Za-z0-9_-]", "_", lipid_name)

        tryCatch(
          {
            # Extract data for this lipid
            if (lipid_name %in% rownames(abundance_data)) {
              lipid_values <- abundance_data[lipid_name, ]

              # Create data frame for plotting
              plot_data <- data.frame(
                sample = colnames(abundance_data),
                abundance = as.numeric(lipid_values),
                group = col_data$group
              )

              # Remove NA values
              plot_data <- plot_data[!is.na(plot_data$abundance), ]

              # Check if we have enough data in each group
              group_counts <- table(plot_data$group)
              if (all(group_counts >= 2)) {
                # Get p-value for this lipid
                pval <- sig_lipids[[pval_col]][i]
                pval_text <- if (pval < 0.001) {
                  "p < 0.001"
                } else {
                  sprintf("p = %.3f", pval)
                }

                # Determine significance stars
                sig_label <- if (pval < 0.001) {
                  "***"
                } else if (pval < 0.01) {
                  "**"
                } else if (pval < 0.05) {
                  "*"
                } else {
                  "ns"
                }

                # Create boxplot with ggplot2
                p <- ggplot(plot_data, aes(x = group, y = abundance, fill = group)) +
                  geom_boxplot(outlier.shape = NA) +
                  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
                  geom_signif(
                    comparisons = list(c("healthy", "severe")),
                    annotations = sig_label,
                    textsize = 6,
                    vjust = 0.5,
                    map_signif_level = FALSE
                  ) +
                  labs(
                    title = lipid_name,
                    subtitle = pval_text,
                    x = "Group",
                    y = "Log10 Abundance"
                  ) +
                  theme_bw() +
                  theme(
                    plot.title = element_text(hjust = 0.5, face = "bold"),
                    plot.subtitle = element_text(hjust = 0.5),
                    legend.position = "none"
                  ) +
                  scale_fill_manual(values = c("healthy" = "#4DAF4A", "severe" = "#E41A1C"))

                # Save plot
                ggsave(
                  filename = file.path(
                    base_dir,
                    "02.DiffExp/03.Individual_lipid_boxplots",
                    paste0(sprintf("%03d", i), "_", safe_name, ".png")
                  ),
                  plot = p,
                  width = 8, height = 6, dpi = 120
                )

                success_count <- success_count + 1

                if (i %% 10 == 0) {
                  cat("    Progress:", i, "/", max_plots, "(", success_count, "successful )\n")
                }
              } else {
                failed_count <- failed_count + 1
                cat("    Warning: Not enough samples for", lipid_name,
                    "(healthy:", group_counts["healthy"], ", severe:", group_counts["severe"], ")\n")
              }
            } else {
              failed_count <- failed_count + 1
              cat("    Warning: Lipid", lipid_name, "not found in data\n")
            }
          },
          error = function(e) {
            failed_count <<- failed_count + 1
            cat(
              "    Warning: Could not create plot for", lipid_name, "-",
              conditionMessage(e), "\n"
            )
          }
        )
      }

      cat("  Boxplots created:", success_count, "/", max_plots, "\n")
      cat("  Failed:", failed_count, "\n")

      cat("  ✓ Individual lipid boxplots complete\n")
    } else {
      cat("  No significant lipids found\n")
    }
  } else {
    cat("  Warning: Could not find p-value column\n")
  }
} else {
  cat("  Warning: No results available for boxplot creation\n")
}


##########################################################################################
# 7 & 8. LIPID SET ENRICHMENT ANALYSIS (LSEA) & OVER REPRESENTATION ANALYSIS (ORA)
##########################################################################################

cat("\n=== 7 & 8. LSEA / ORA ===\n")

run_lsea_ora(
  deSp_se,
  char_list = list(by_class = "class", by_category = "Category"),
  lsea_dir = file.path(base_dir, "03.Enrichment/LSEA"),
  ora_dir = file.path(base_dir, "03.Enrichment/ORA")
)


##########################################################################################
# ANALYSIS COMPLETE
##########################################################################################

cat("\n========================================\n")
cat("POST-HOC ANALYSIS HEALTHY vs SEVERE COMPLETED!\n")
cat("========================================\n")
