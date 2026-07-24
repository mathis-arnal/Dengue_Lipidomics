##########################################################################################
# LIPIDOMICS ANALYSIS: THREE GROUPS (HEALTHY + MILD + SEVERE)
##########################################################################################
#
# COMPLETE ANALYSIS WORKFLOW:
# 1. Data Loading & Preparation
# 2. Profiling (Data Quality, Cross-Sample Variability, Dimensionality Reduction)
# 3. Correlation Analysis
# 4. Lipid Characteristics
# 5. Differential Expression Analysis (ANOVA)
# 6. Lipid Set Enrichment Analysis (LSEA)
# 7. Over Representation Analysis (ORA)
#
# This script uses ALL capabilities of LipidSigR for comprehensive 3-group analysis
##########################################################################################

# Load required libraries
library(LipidSigR)
library(SummarizedExperiment)
library(dplyr)
library(readr)
library(rgoslin)
library(ggplot2)

source("code/lipidomics_functions.R")

# Set working directory
setwd("E:/Dengue_lipidomics")

# Define base directory for outputs
base_dir <- "E:/Dengue_lipidomics/analysis/PBMC/Three_groups/D0"

cat("\n========================================\n")
cat("THREE-GROUP LIPIDOMICS ANALYSIS\n")
cat("Groups: Healthy, Mild, Severe\n")
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
  "02.DiffExp/01.ANOVA",
  "02.DiffExp/02.Visualizations",
  "02.DiffExp/03.Individual_lipid_boxplots",
  "03.Enrichment/LSEA",
  "03.Enrichment/ORA"
))

cat("Loading group information and abundance data...\n")
loaded <- load_and_match_samples(
  "data/PBMC/group_information_table_healthy_vs_sick_patients_D0.tsv",
  "data/PBMC/PBMC_lipid_abundance_data_D0.tsv"
)
abundance <- loaded$abundance
group_info <- loaded$group_info
cat(" The 6 missing samples are expected, KT-312D0 doesn't have enough PBMC.",
    "Les 5 healthy ne rentraient pas dans la machine.\n")

annotated <- clean_and_annotate_lipids(abundance)

# Construct SummarizedExperiment object + process
se_result <- build_processed_se(
  annotated$abundance_filtered, annotated$goslin_annotation, group_info,
  se_type = "de_multiple", exclude_missing_pct = 70, paired_sample = NULL
)
se <- se_result$se
processed_se <- se_result$processed_se


##########################################################################################
# 2. PROFILING - DATA QUALITY
##########################################################################################

cat("\n=== 2. PROFILING - DATA QUALITY ===\n")

# Plot data processing effects
cat("Generating data quality plots...\n")
data_process_plots <- plot_data_process(se, processed_se)

save_static_png(data_process_plots$static_boxPlot_before,
  file.path(base_dir, "01.Profiling/00.Data_quality/BoxPlot_before_process.png"),
  width = 1400, height = 1000
)
save_static_png(data_process_plots$static_boxPlot_after,
  file.path(base_dir, "01.Profiling/00.Data_quality/BoxPlot_after_process.png"),
  width = 1400, height = 1000
)
save_static_png(data_process_plots$static_densityPlot_before,
  file.path(base_dir, "01.Profiling/00.Data_quality/DensityPlot_before_process.png"),
  width = 1400, height = 1000
)
save_static_png(data_process_plots$static_densityPlot_after,
  file.path(base_dir, "01.Profiling/00.Data_quality/DensityPlot_after_process.png"),
  width = 1400, height = 1000
)

cat("  ✓ Data quality plots saved\n")


# Cross-sample variability
cat("Analyzing cross-sample variability...\n")
cross_var_result <- cross_sample_variability(se)

save_static_png(cross_var_result$static_lipid_number_barPlot,
  file.path(base_dir, "01.Profiling/01.Cross-sample_variability/Expressed_lipid_numbers.png"),
  width = 1600, height = 1000
)
save_static_png(cross_var_result$static_lipid_amount_barPlot,
  file.path(base_dir, "01.Profiling/01.Cross-sample_variability/Histogram_lipid_amount.png"),
  width = 1400, height = 1000
)
save_static_png(cross_var_result$static_lipid_distribution,
  file.path(base_dir, "01.Profiling/01.Cross-sample_variability/Lipid_abundance_distribution.png"),
  width = 1400, height = 1000
)

cat("  ✓ Cross-sample variability analysis complete\n")


##########################################################################################
# 3. PROFILING - DIMENSIONALITY REDUCTION
##########################################################################################

cat("\n=== 3. PROFILING - DIMENSIONALITY REDUCTION ===\n")

# PCA Analysis
cat("Running PCA analysis...\n")
result_pca <- dr_pca(
  processed_se,
  scaling = TRUE,
  centering = TRUE,
  clustering = "kmeans",
  cluster_num = 3,
  kmedoids_metric = NULL,
  distfun = NULL,
  hclustfun = NULL,
  eps = NULL,
  minPts = NULL,
  feature_contrib_pc = c(1, 2),
  plot_topN = 10
)

# Severity-colored PCA (color = true clinical group, shape = unsupervised k-means
# cluster, 95% ellipse per group) -- see plot_pca_by_severity() in
# lipidomics_functions.R for why clustering="group_info" can't be used directly here.
pca_severity <- plot_pca_by_severity(result_pca, group_info)

save_static_png(pca_severity$plot,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/PCA_plot_clinical_group.png"),
  width = 1400, height = 1000
)
save_interactive_rds(pca_severity$plot_interactive,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/PCA_plot_clinical_group_interactive.rds"))

save_static_png(result_pca$static_pca,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/PCA_plot.png"),
  width = 1400, height = 1000
)
save_interactive_rds(result_pca$interactive_pca,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/PCA_plot_interactive.rds"))

save_static_png(result_pca$static_screePlot,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Scree_plot.png"),
  width = 1400, height = 1000
)
save_interactive_rds(result_pca$interactive_screePlot,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Scree_plot_interactive.rds"))

save_static_png(result_pca$static_feature_contribution,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Feature_contribution.png"),
  width = 1400, height = 1000
)
save_interactive_rds(result_pca$interactive_feature_contribution,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Feature_contribution_interactive.rds"))

save_static_png(result_pca$static_variablePlot,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Variable_correlation.png"),
  width = 1400, height = 1000
)
save_interactive_rds(result_pca$interactive_variablePlot,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Variable_correlation_interactive.rds"))

write_tsv(result_pca$table_pca_contribution,
  file = file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/table_pca_contribution.tsv")
)
write_tsv(result_pca$pca_rotated_data,
  file = file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/table_pca_rotated_data.tsv")
)

cat("  ✓ PCA complete\n")


# t-SNE Analysis
cat("Running t-SNE analysis...\n")
result_tsne <- dr_tsne(
  processed_se,
  pca = TRUE,
  perplexity = 5,
  max_iter = 500,
  clustering = "kmeans",
  cluster_num = 3,
  kmedoids_metric = NULL,
  distfun = NULL,
  hclustfun = NULL,
  eps = NULL,
  minPts = NULL
)

save_static_png(result_tsne$static_tsne,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/t-SNE/t-SNE_plot.png"),
  width = 1400, height = 1000
)
save_interactive_rds(result_tsne$interactive_tsne,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/t-SNE/t-SNE_plot_interactive.rds"))

write_tsv(result_tsne$tsne_result,
  file = file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/t-SNE/table_tsne_results.tsv")
)

cat("  ✓ t-SNE complete\n")


# UMAP Analysis
cat("Running UMAP analysis...\n")
result_umap <- dr_umap(
  processed_se,
  n_neighbors = 15,
  scaling = TRUE,
  umap_metric = "euclidean",
  clustering = "kmeans",
  cluster_num = 3,
  kmedoids_metric = NULL,
  distfun = NULL,
  hclustfun = NULL,
  eps = NULL,
  minPts = NULL
)

save_static_png(result_umap$static_umap,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/UMAP/UMAP_plot.png"),
  width = 1400, height = 1000
)
save_interactive_rds(result_umap$interactive_umap,
  file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/UMAP/UMAP_plot_interactive.rds"))

write_tsv(result_umap$umap_result,
  file = file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/UMAP/UMAP_results.tsv")
)

cat("  ✓ UMAP complete\n")


##########################################################################################
# 4. CORRELATION ANALYSIS
##########################################################################################

cat("\n=== 4. CORRELATION ANALYSIS ===\n")

# Correlation by sample
cat("Computing correlation heatmap by sample...\n")
result_heatmap_sample <- heatmap_correlation(
  processed_se,
  char = NULL,
  transform = "log10",
  correlation = "spearman",
  distfun = "euclidean",
  hclustfun = "ward.D2",
  type = "sample"
)

save_static_png(result_heatmap_sample$static_heatmap,
  file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_samples/Heatmap_by_samples.png"),
  width = 2000, height = 2000
)
save_interactive_rds(result_heatmap_sample$interactive_heatmap,
  file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_samples/Heatmap_by_samples_interactive.rds"))

cat("  ✓ Sample correlation complete\n")


# Correlation by category
cat("Computing correlation heatmap by category...\n")
result_heatmap_category <- heatmap_correlation(
  processed_se,
  char = "Category",
  transform = "log10",
  correlation = "spearman",
  distfun = "euclidean",
  hclustfun = "ward.D2",
  type = "class"
)

save_static_png(result_heatmap_category$static_heatmap,
  file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_category/Heatmap_by_category.png"),
  width = 1500, height = 1500
)
save_interactive_rds(result_heatmap_category$interactive_heatmap,
  file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_category/Heatmap_by_category_interactive.rds"))

cat("  ✓ Category correlation complete\n")


# Correlation by lipid class
cat("Computing correlation heatmap by lipid class...\n")
result_heatmap_class <- heatmap_correlation(
  processed_se,
  char = "class",
  transform = "log10",
  correlation = "spearman",
  distfun = "euclidean",
  hclustfun = "ward.D2",
  type = "class"
)

save_static_png(result_heatmap_class$static_heatmap,
  file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_class/Heatmap_by_class.png"),
  width = 1500, height = 1500
)
save_interactive_rds(result_heatmap_class$interactive_heatmap,
  file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_class/Heatmap_by_class_interactive.rds"))

write_tsv(as.data.frame(result_heatmap_class$corr_coef_matrix),
  file = file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_class/Correlation_matrix.tsv")
)

cat("  ✓ Class correlation complete\n")


##########################################################################################
# 5. LIPID CHARACTERISTICS PER SAMPLE
##########################################################################################

cat("\n=== 5. LIPID CHARACTERISTICS ===\n")

# Lipid profiling by class
cat("Analyzing lipid composition by class...\n")
result_lipid_class <- lipid_profiling(processed_se, char = "class")

save_static_png(result_lipid_class$static_char_barPlot,
  file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Lipid_class_per_sample.png"),
  width = 1000, height = 1000
)
save_static_png(result_lipid_class$static_lipid_composition,
  file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Lipid_composition_per_sample.png"),
  width = 1000, height = 1000
)


# Lipid profiling by Total Carbon
cat("Analyzing lipid composition by Total.C...\n")
result_totalC <- lipid_profiling(processed_se, char = "Total.C")

save_static_png(result_totalC$static_char_barPlot,
  file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Total_C_per_sample.png"),
  width = 1000, height = 1000
)


# Lipid profiling by Total Double Bonds
cat("Analyzing lipid composition by Total.DB...\n")
result_totalDB <- lipid_profiling(processed_se, char = "Total.DB")

save_static_png(result_totalDB$static_char_barPlot,
  file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Total_DB_per_sample.png"),
  width = 1000, height = 1000
)


# Lipid profiling by Category
cat("Analyzing lipid composition by Category...\n")
result_category <- lipid_profiling(processed_se, char = "Category")

save_static_png(result_category$static_char_barPlot,
  file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Category_per_sample.png"),
  width = 1000, height = 1000
)

cat("  ✓ Lipid characteristics analysis complete\n")


##########################################################################################
# 6. DIFFERENTIAL EXPRESSION ANALYSIS - ANOVA
##########################################################################################

cat("\n=== 6. DIFFERENTIAL EXPRESSION ANALYSIS (ANOVA) ===\n")

# Multi-group differential expression using ANOVA
cat("Running ANOVA for 3-group comparison...\n")
deSp_se <- deSp_multiGroup(
  processed_se,
  test = "One-way ANOVA",
  ref_group = "healthy",
  significant = "pval",
  p_cutoff = 0.05,
  transform = "log10"
)

# Extract results
res_list <- extract_summarized_experiment(deSp_se)

# Check the structure of results
cat("  Result structure:\n")
print(names(res_list))

# Get the differential expression results
# For multi-group ANOVA, results are in all_deSp_result and sig_deSp_result
if (!is.null(res_list$all_deSp_result)) {
  desp_df <- as.data.frame(res_list$all_deSp_result)

  # Count significant lipids (deSp_multiGroup() output uses the same lowercase
  # pval/padj column names as deSp_twoGroup(), not limma-style P.Value/adj.P.Val)
  n_sig <- sum(desp_df$padj < 0.05, na.rm = TRUE)
  cat("  Total lipids analyzed:", nrow(desp_df), "\n")
  cat("  Significant lipids (padj < 0.05):", n_sig, "\n")

  # Save all results
  write_tsv(desp_df,
    file = file.path(base_dir, "02.DiffExp/01.ANOVA/ANOVA_all_results.tsv")
  )

  # Save significant results separately
  if (n_sig > 0) {
    sig_df <- desp_df %>% filter(padj < 0.05)
    write_tsv(sig_df,
      file = file.path(base_dir, "02.DiffExp/01.ANOVA/ANOVA_significant_results.tsv")
    )
  }

  cat("  Results extracted and saved.\n")
} else {
  cat("  Warning: Could not find expected result structure.\n")
  cat("  Available result fields:\n")
  print(names(res_list))
  desp_df <- NULL
}

# Plot differential expression results
deSp_plot <- plot_deSp_multiGroup(deSp_se)

# Lollipop chart of significant lipids
save_static_png(deSp_plot$static_de_lipid,
  file.path(base_dir, "02.DiffExp/02.Visualizations/Lollipop_chart.png"),
  width = 1000, height = 1000
)

cat("  ✓ ANOVA analysis complete\n")


# Pairwise comparisons within ANOVA (if available in LipidSigR)
# This provides detailed comparisons between each pair of groups
cat("Generating pairwise comparison plots...\n")

# Generate heatmap of significantly different lipids
save_static_png(deSp_plot$static_heatmap,
  file.path(base_dir, "02.DiffExp/02.Visualizations/Heatmap_significant_lipids.png"),
  width = 1500, height = 1200
)
save_interactive_rds(deSp_plot$interactive_heatmap,
  file.path(base_dir, "02.DiffExp/02.Visualizations/Heatmap_significant_lipids_interactive.rds"))

cat("  ✓ Visualization complete\n")


# Individual lipid species boxplots
cat("\n--- Creating individual lipid abundance boxplots ---\n")

# Get significant lipids from ANOVA results
if (exists("desp_df") && !is.null(desp_df) && nrow(desp_df) > 0) {
  # Check which p-value column is available
  # Nominal (uncorrected) p-value first -- boxplots should cover every nominally
  # significant lipid, not just the (usually much smaller/empty) FDR-significant set.
  # deSp_multiGroup() output uses the same lowercase pval/padj columns as
  # deSp_twoGroup(), not limma-style P.Value/adj.P.Val -- this previously meant the
  # column lookup below always failed and no three-group boxplots were ever generated.
  pval_col <- if ("pval" %in% colnames(desp_df)) {
    "pval"
  } else if ("padj" %in% colnames(desp_df)) {
    "padj"
  } else {
    NULL
  }

  if (!is.null(pval_col)) {
    # Filter significant lipids (p < 0.05)
    sig_lipids <- desp_df %>%
      filter(.data[[pval_col]] < 0.05) %>%
      arrange(.data[[pval_col]])

    n_sig <- nrow(sig_lipids)
    cat("  Creating boxplots for", n_sig, "significant lipids...\n")

    if (n_sig > 0) {
      # Limit to top 50 most significant to avoid too many plots
      max_plots <- min(50, n_sig)
      cat("  (Showing top", max_plots, "most significant lipids)\n")

      for (i in 1:max_plots) {
        # Get lipid name - check if column exists
        lipid_name <- if ("feature" %in% colnames(sig_lipids)) {
          sig_lipids$feature[i]
        } else if ("lipid" %in% colnames(sig_lipids)) {
          sig_lipids$lipid[i]
        } else {
          rownames(sig_lipids)[i]
        }

        # Create safe filename (replace special characters)
        safe_name <- gsub("[^A-Za-z0-9_-]", "_", lipid_name)

        tryCatch({
          # Create boxplot for this lipid across all three groups
          boxplot_result <- boxPlot_feature_multiGroup(
            processed_se,
            feature = lipid_name,
            test = "One-way ANOVA",
            transform = "log10"
          )

          save_static_png(boxplot_result$static_boxPlot,
            file.path(
              base_dir,
              "02.DiffExp/03.Individual_lipid_boxplots",
              paste0(sprintf("%03d", i), "_", safe_name, ".png")
            ),
            width = 800, height = 600, res = 120
          )

          # Progress indicator every 10 lipids
          if (i %% 10 == 0) {
            cat("    Progress:", i, "/", max_plots, "\n")
          }
        }, error = function(e) {
          cat("    Warning: Could not create plot for", lipid_name, "-",
              conditionMessage(e), "\n")
        })
      }

      cat("  ✓ Individual lipid boxplots complete\n")
    } else {
      cat("  No significant lipids found (p < 0.05)\n")
    }
  } else {
    cat("  Warning: Could not find p-value column in results\n")
    cat("  Available columns:", paste(colnames(desp_df), collapse = ", "), "\n")
  }
} else {
  cat("  Warning: ANOVA results not available for boxplot creation\n")
}


##########################################################################################
# 7 & 8. LIPID SET ENRICHMENT ANALYSIS (LSEA) & OVER REPRESENTATION ANALYSIS (ORA)
##########################################################################################

cat("\n=== 7 & 8. LSEA / ORA ===\n")

run_lsea_ora(
  deSp_se,
  char_list = list(all_characteristics = NULL, by_class = "class", by_category = "Category"),
  lsea_dir = file.path(base_dir, "03.Enrichment/LSEA"),
  ora_dir = file.path(base_dir, "03.Enrichment/ORA")
)


##########################################################################################
# ANALYSIS COMPLETE
##########################################################################################

cat("\n========================================\n")
cat("THREE-GROUP ANALYSIS COMPLETED!\n")
cat("========================================\n")
cat("\nOutput directory:", base_dir, "\n")
cat("\nGenerated outputs:\n")
cat("  - Profiling plots (data quality, PCA, t-SNE, UMAP)\n")
cat("  - Correlation heatmaps\n")
cat("  - Lipid characteristics profiles\n")
cat("  - ANOVA differential expression results\n")
cat("  - LSEA enrichment analysis\n")
cat("\n")
