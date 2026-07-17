##########################################################################################
# LIPIDOMICS ANALYSIS: POST-HOC MILD vs SEVERE
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
# This script performs pairwise comparison: Mild vs Severe
##########################################################################################

# Load required libraries
library(LipidSigR)
library(SummarizedExperiment)
library(dplyr)
library(readr)
library(rgoslin)
library(ggplot2)
library(ggsignif)

# Set working directory
setwd("/Users/loictalignani/research/project/lipidomics")

# Define base directory for outputs
base_dir <- "/Users/loictalignani/research/project/lipidomics/analysis/PostHoc_Mild_vs_Severe/D0"

cat("\n========================================\n")
cat("POST-HOC ANALYSIS: MILD vs SEVERE\n")
cat("========================================\n")


##########################################################################################
# 1. DATA LOADING & PREPARATION
##########################################################################################

cat("\n=== 1. DATA LOADING & PREPARATION ===\n")

# Create output directories
dir.create(file.path(base_dir, "01.Profiling/00.Data_quality"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "01.Profiling/01.Cross-sample_variability"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/t-SNE"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/UMAP"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_samples"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_category"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_class"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "01.Profiling/04.Lipid_characteristics"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "02.DiffExp/01.t-test_results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "02.DiffExp/02.Visualizations/Volcano"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "02.DiffExp/02.Visualizations/MA_plot"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "02.DiffExp/02.Visualizations/Heatmap"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "02.DiffExp/03.Individual_lipid_boxplots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "03.Enrichment/LSEA"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "03.Enrichment/ORA"), recursive = TRUE, showWarnings = FALSE)

# Load group information
cat("Loading group information...\n")
group_info <- read_tsv("data/lipidsig_datasets/healthy_vs_sick_patients/group_information_table_MILD_vs_SEVERE_patients_D0.tsv",
  show_col_types = FALSE
)
cat("  Groups found:", unique(group_info$group), "\n")
cat("  Sample count:", nrow(group_info), "\n")

# Load abundance data
cat("Loading abundance data...\n")
abundance <- read_tsv("data/lipidsig_datasets/healthy_vs_sick_patients/MILD_SEVERE_lipidomics.tsv",
  show_col_types = FALSE
)
cat("  Features:", nrow(abundance), "\n")
cat("  Samples:", ncol(abundance) - 1, "\n")

# Parse lipid names with rgoslin
cat("Parsing lipid nomenclature with rgoslin...\n")
parse_lipid <- rgoslin::parseLipidNames(lipidNames = abundance$feature)

# Filter recognized lipids
recognized_lipids <- parse_lipid$Original.Name[which(parse_lipid$Grammar != "NOT_PARSEABLE")]
abundance_filtered <- abundance %>% dplyr::filter(feature %in% recognized_lipids)
goslin_annotation <- parse_lipid %>% dplyr::filter(Original.Name %in% recognized_lipids)

cat("  Recognized lipids:", length(recognized_lipids), "/", nrow(abundance), "\n")
cat("  Recognition rate:", round(100 * length(recognized_lipids) / nrow(abundance), 1), "%\n")

# Construct SummarizedExperiment object
cat("Creating SummarizedExperiment object...\n")
se <- as_summarized_experiment(
  abundance_filtered,
  goslin_annotation,
  group_info = group_info,
  se_type = "de_two", # Two-group comparison
  paired_sample = FALSE
)

cat("  ✓ SummarizedExperiment created\n")

# Data processing
cat("Processing data (filtering, normalization, transformation)...\n")
processed_se <- data_process(
  se,
  exclude_missing = TRUE,
  exclude_missing_pct = 50,
  replace_na_method = "min",
  replace_na_method_ref = 0.5,
  normalization = "Percentage",
  transform = "log10"
)

cat("  ✓ Data processing complete\n")


##########################################################################################
# 2. PROFILING
##########################################################################################

cat("\n=== 2. PROFILING ===\n")

# Data quality assessment - compare before/after processing
cat("Assessing data quality...\n")
plot_quality <- plot_data_process(se, processed_se)

png(file.path(base_dir, "01.Profiling/00.Data_quality/Data_quality_boxplot_before.png"),
  width = 1200, height = 800, res = 150
)
print(plot_quality$static_boxPlot_before)
dev.off()

png(file.path(base_dir, "01.Profiling/00.Data_quality/Data_quality_boxplot_after.png"),
  width = 1200, height = 800, res = 150
)
print(plot_quality$static_boxPlot_after)
dev.off()

png(file.path(base_dir, "01.Profiling/00.Data_quality/Data_quality_density_before.png"),
  width = 1200, height = 800, res = 150
)
print(plot_quality$static_densityPlot_before)
dev.off()

png(file.path(base_dir, "01.Profiling/00.Data_quality/Data_quality_density_after.png"),
  width = 1200, height = 800, res = 150
)
print(plot_quality$static_densityPlot_after)
dev.off()

# Cross-sample variability
cat("Analyzing cross-sample variability...\n")
plot_variability <- cross_sample_variability(processed_se)

png(file.path(base_dir, "01.Profiling/01.Cross-sample_variability/Expressed_lipid_numbers.png"),
  width = 1200, height = 800, res = 150
)
print(plot_variability$static_lipid_number_barPlot)
dev.off()

png(file.path(base_dir, "01.Profiling/01.Cross-sample_variability/Lipid_abundance_distribution.png"),
  width = 1200, height = 800, res = 150
)
print(plot_variability$static_lipid_distribution)
dev.off()

# PCA
cat("Performing PCA...\n")
plot_pca <- dr_pca(
  processed_se,
  scaling = TRUE,
  centering = TRUE,
  clustering = "kmeans",
  cluster_num = 2
)

png(file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/PCA_plot.png"),
  width = 1000, height = 800, res = 150
)
print(plot_pca$static_pca)
dev.off()

png(file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Scree_plot.png"),
  width = 1000, height = 800, res = 150
)
print(plot_pca$static_screePlot)
dev.off()

# Save interactive PCA plots as RDS
if (!is.null(plot_pca$interactive_pca)) {
  saveRDS(
    plot_pca$interactive_pca,
    file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/PCA_plot_interactive.rds")
  )
}

if (!is.null(plot_pca$interactive_screePlot)) {
  saveRDS(
    plot_pca$interactive_screePlot,
    file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Scree_plot_interactive.rds")
  )
}

if (!is.null(plot_pca$interactive_variablePlot)) {
  saveRDS(
    plot_pca$interactive_variablePlot,
    file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Variable_correlation_interactive.rds")
  )
}

# t-SNE
cat("Performing t-SNE...\n")
plot_tsne <- dr_tsne(
  processed_se,
  pca = TRUE,
  perplexity = 5,
  clustering = "kmeans",
  cluster_num = 2
)

png(file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/t-SNE/tSNE_plot.png"),
  width = 1000, height = 800, res = 150
)
print(plot_tsne$static_tsne)
dev.off()

# Save interactive t-SNE plot as RDS
if (!is.null(plot_tsne$interactive_tsne)) {
  saveRDS(
    plot_tsne$interactive_tsne,
    file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/t-SNE/t-SNE_plot_interactive.rds")
  )
}

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

png(file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/UMAP/UMAP_plot.png"),
  width = 1000, height = 800, res = 150
)
print(plot_umap$static_umap)
dev.off()

# Save interactive UMAP plot as RDS
if (!is.null(plot_umap$interactive_umap)) {
  saveRDS(
    plot_umap$interactive_umap,
    file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/UMAP/UMAP_plot_interactive.rds")
  )
}

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

png(file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_samples/Correlation_heatmap_samples.png"),
  width = 1200, height = 1000, res = 150
)
print(plot_corr_samples$static_heatmap)

# Save interactive heatmap as RDS
if (!is.null(plot_corr_samples$interactive_heatmap)) {
  saveRDS(
    plot_corr_samples$interactive_heatmap,
    file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_samples/Correlation_heatmap_samples_interactive.rds")
  )
}
dev.off()

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

png(file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_category/Correlation_heatmap_category.png"),
  width = 1200, height = 1000, res = 150
)
print(plot_corr_category$static_heatmap)

# Save interactive heatmap as RDS
if (!is.null(plot_corr_category$interactive_heatmap)) {
  saveRDS(
    plot_corr_category$interactive_heatmap,
    file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_category/Correlation_heatmap_category_interactive.rds")
  )
}
dev.off()

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

png(file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_class/Correlation_heatmap_class.png"),
  width = 1500, height = 1200, res = 150
)
print(plot_corr_class$static_heatmap)

# Save interactive heatmap as RDS
if (!is.null(plot_corr_class$interactive_heatmap)) {
  saveRDS(
    plot_corr_class$interactive_heatmap,
    file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_class/Correlation_heatmap_class_interactive.rds")
  )
}
dev.off()

cat("  ✓ Correlation analysis complete\n")


##########################################################################################
# 4. LIPID CHARACTERISTICS
##########################################################################################

cat("\n=== 4. LIPID CHARACTERISTICS ===\n")

# Characteristics by class
cat("Analyzing lipid characteristics by class...\n")
result_class <- lipid_profiling(processed_se, char = "class")

png(file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Lipid_characteristics_by_class.png"),
  width = 1200, height = 1000, res = 150
)
print(result_class$static_char_barPlot)
dev.off()

png(file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Lipid_composition_per_sample.png"),
  width = 1000, height = 1000, res = 150
)
print(result_class$static_lipid_composition)
dev.off()

# Characteristics by Total.C
cat("Analyzing lipid characteristics by Total.C...\n")
result_totalC <- lipid_profiling(processed_se, char = "Total.C")

png(file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Lipid_characteristics_by_TotalC.png"),
  width = 1200, height = 1000, res = 150
)
print(result_totalC$static_char_barPlot)
dev.off()

# Characteristics by Total.DB
cat("Analyzing lipid characteristics by Total.DB...\n")
result_totalDB <- lipid_profiling(processed_se, char = "Total.DB")

png(file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Lipid_characteristics_by_TotalDB.png"),
  width = 1200, height = 1000, res = 150
)
print(result_totalDB$static_char_barPlot)
dev.off()

# Characteristics by Category
cat("Analyzing lipid characteristics by Category...\n")
result_category <- lipid_profiling(processed_se, char = "Category")

png(file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Lipid_characteristics_by_Category.png"),
  width = 1200, height = 1000, res = 150
)
print(result_category$static_char_barPlot)
dev.off()

cat("  ✓ Lipid characteristics analysis complete\n")


##########################################################################################
# 5. DIFFERENTIAL EXPRESSION ANALYSIS (T-TEST)
##########################################################################################

cat("\n=== 5. DIFFERENTIAL EXPRESSION ANALYSIS (t-test) ===\n")

# Two-group differential expression using t-test
cat("Running t-test: Mild vs Severe...\n")
deSp_se <- deSp_twoGroup(
  processed_se,
  ref_group = "mild",
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

# Volcano plot - Static version
png(file.path(base_dir, "02.DiffExp/02.Visualizations/Volcano/Volcano_plot.png"),
  width = 1200, height = 1000, res = 150
)
print(deSp_plot$static_volcanoPlot)
dev.off()

# Volcano plot - Interactive version
if (!is.null(deSp_plot$interactive_volcanoPlot)) {
  htmlwidgets::saveWidget(
    deSp_plot$interactive_volcanoPlot,
    file.path(base_dir, "02.DiffExp/02.Visualizations/Volcano/Volcano_plot_interactive.html"),
    selfcontained = TRUE
  )
  # Also save as RDS for embedding in Rmd
  saveRDS(
    deSp_plot$interactive_volcanoPlot,
    file.path(base_dir, "02.DiffExp/02.Visualizations/Volcano/Volcano_plot_interactive.rds")
  )
  cat("  Interactive volcano plot saved\n")
}

# MA plot - Static version
png(file.path(base_dir, "02.DiffExp/02.Visualizations/MA_plot/MA_plot.png"),
  width = 1200, height = 1000, res = 150
)
print(deSp_plot$static_maPlot)
dev.off()

# MA plot - Interactive version
if (!is.null(deSp_plot$interactive_maPlot)) {
  htmlwidgets::saveWidget(
    deSp_plot$interactive_maPlot,
    file.path(base_dir, "02.DiffExp/02.Visualizations/MA_plot/MA_plot_interactive.html"),
    selfcontained = TRUE
  )
  # Also save as RDS for embedding in Rmd
  saveRDS(
    deSp_plot$interactive_maPlot,
    file.path(base_dir, "02.DiffExp/02.Visualizations/MA_plot/MA_plot_interactive.rds")
  )
  cat("  Interactive MA plot saved\n")
}

# Lollipop chart
png(file.path(base_dir, "02.DiffExp/02.Visualizations/Lollipop_chart.png"),
  width = 1200, height = 1000, res = 150
)
print(deSp_plot$static_de_lipid)
dev.off()

# Heatmap
if (!is.null(deSp_plot$static_heatmap)) {
  png(file.path(base_dir, "02.DiffExp/02.Visualizations/Heatmap/Heatmap_significant_lipids.png"),
    width = 1500, height = 1200, res = 150
  )
  print(deSp_plot$static_heatmap)
  dev.off()
}

cat("  ✓ Differential expression analysis complete\n")


##########################################################################################
# 5A. PLS-DA ANALYSIS
##########################################################################################

cat("\n=== 5A. PLS-DA ANALYSIS ===\n")

# Create directory for PLS-DA outputs
dir.create(file.path(base_dir, "02.DiffExp/02.Visualizations/PLS-DA"), recursive = TRUE, showWarnings = FALSE)

# Perform PLS-DA
cat("Running PLS-DA...\n")
result_plsda <- dr_plsda(
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
)

# Save PLS-DA results tables
write.csv(result_plsda$plsda_result,
  file = file.path(base_dir, "02.DiffExp/02.Visualizations/PLS-DA/Table_of_sample_variate.csv"),
  row.names = FALSE
)

write.csv(result_plsda$table_plsda_loading,
  file = file.path(base_dir, "02.DiffExp/02.Visualizations/PLS-DA/Table_of_sample_loading.csv"),
  row.names = FALSE
)

# Save static PLS-DA plot
png(file.path(base_dir, "02.DiffExp/02.Visualizations/PLS-DA/PLSDA_plot.png"),
  width = 1200, height = 800, res = 150
)
print(result_plsda$static_plsda)
dev.off()

# Save interactive PLS-DA plot
if (!is.null(result_plsda$interacitve_plsda)) {
  saveRDS(
    result_plsda$interacitve_plsda,
    file.path(base_dir, "02.DiffExp/02.Visualizations/PLS-DA/PLSDA_plot_interactive.rds")
  )
}

# Save static loading plot
png(file.path(base_dir, "02.DiffExp/02.Visualizations/PLS-DA/Loading_plot.png"),
  width = 800, height = 800, res = 150
)
print(result_plsda$static_loadingPlot)
dev.off()

# Save interactive loading plot
if (!is.null(result_plsda$interactive_loadingPlot)) {
  saveRDS(
    result_plsda$interactive_loadingPlot,
    file.path(base_dir, "02.DiffExp/02.Visualizations/PLS-DA/Loading_plot_interactive.rds")
  )
}

cat("  ✓ PLS-DA analysis complete\n")


##########################################################################################
# 5B. DIFFERENTIAL EXPRESSION OF LIPID CHARACTERISTICS
##########################################################################################

cat("\n=== 5B. DIFFERENTIAL EXPRESSION OF LIPID CHARACTERISTICS ===\n")

# Analyze differential expression of lipid characteristics
cat("Analyzing differential expression of lipid characteristics...\n")

# By Total Carbon (Total.C)
cat("  - By Total Carbon\n")
deChar_totalC_se <- deChar_twoGroup(
  processed_se,
  char = "Total.C",
  ref_group = "mild",
  test = "t-test",
  significant = "pval",
  p_cutoff = 0.05,
  FC_cutoff = 1,
  transform = "log10"
)

deChar_totalC_plot <- plot_deChar_twoGroup(deChar_totalC_se)

png(file.path(base_dir, "02.DiffExp/02.Visualizations/deChar_TotalC_barplot.png"),
  width = 1200, height = 1000, res = 150
)
print(deChar_totalC_plot$static_char_barPlot)
dev.off()

# By Category
cat("  - By Category\n")
deChar_category_se <- deChar_twoGroup(
  processed_se,
  char = "Category",
  ref_group = "mild",
  test = "t-test",
  significant = "pval",
  p_cutoff = 0.05,
  FC_cutoff = 1,
  transform = "log10"
)

deChar_category_plot <- plot_deChar_twoGroup(deChar_category_se)

png(file.path(base_dir, "02.DiffExp/02.Visualizations/deChar_Category_barplot.png"),
  width = 1200, height = 1000, res = 150
)
print(deChar_category_plot$static_char_barPlot)
dev.off()

# By Class
cat("  - By Class\n")
deChar_class_se <- deChar_twoGroup(
  processed_se,
  char = "class",
  ref_group = "mild",
  test = "t-test",
  significant = "pval",
  p_cutoff = 0.05,
  FC_cutoff = 1,
  transform = "log10"
)

deChar_class_plot <- plot_deChar_twoGroup(deChar_class_se)

png(file.path(base_dir, "02.DiffExp/02.Visualizations/deChar_Class_barplot.png"),
  width = 1200, height = 1000, res = 150
)
print(deChar_class_plot$static_char_barPlot)
dev.off()

cat("  ✓ Differential expression of lipid characteristics complete\n")


##########################################################################################
# 6. INDIVIDUAL LIPID BOXPLOTS
##########################################################################################

cat("\n=== 6. INDIVIDUAL LIPID BOXPLOTS ===\n")

# Get significant lipids
if (exists("desp_df") && !is.null(desp_df) && nrow(desp_df) > 0) {
  # Check which p-value column is available
  pval_col <- if ("padj" %in% colnames(desp_df)) {
    "padj"
  } else if ("pval" %in% colnames(desp_df)) {
    "pval"
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
                    comparisons = list(c("mild", "severe")),
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
                  scale_fill_manual(values = c("mild" = "#FF7F00", "severe" = "#E41A1C"))

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
                cat(
                  "    Warning: Not enough samples for", lipid_name,
                  "(mild:", group_counts["mild"], ", severe:", group_counts["severe"], ")\n"
                )
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
# 7. LIPID SET ENRICHMENT ANALYSIS (LSEA)
##########################################################################################

cat("\n=== 7. LIPID SET ENRICHMENT ANALYSIS (LSEA) ===\n")

# LSEA by lipid class
cat("Running LSEA by lipid class...\n")
lsea_class_result <- enrichment_lsea(
  deSp_se,
  char = "class",
  rank_by = "statistic",
  significant = "pval",
  p_cutoff = 0.05
)

cat("  LSEA summary (by class):\n")
print(summary(lsea_class_result))

png(file.path(base_dir, "03.Enrichment/LSEA/LSEA_by_class_barplot.png"),
  width = 1200, height = 1000, res = 150
)
print(lsea_class_result$static_barPlot)
dev.off()

cat("  ✓ LSEA by class complete\n\n")

# LSEA by lipid category
cat("Running LSEA by category...\n")
lsea_category_result <- enrichment_lsea(
  deSp_se,
  char = "Category",
  rank_by = "statistic",
  significant = "pval",
  p_cutoff = 0.05
)

cat("  LSEA summary (by category):\n")
print(summary(lsea_category_result))

png(file.path(base_dir, "03.Enrichment/LSEA/LSEA_by_category_barplot.png"),
  width = 1200, height = 1000, res = 150
)
print(lsea_category_result$static_barPlot)
dev.off()

cat("  ✓ LSEA by category complete\n")


##########################################################################################
# 8. OVER REPRESENTATION ANALYSIS (ORA)
##########################################################################################

cat("\n=== 8. OVER REPRESENTATION ANALYSIS (ORA) ===\n")

# ORA by lipid class
cat("Running ORA by lipid class...\n")
ora_class_result <- enrichment_ora(
  deSp_se,
  char = "class",
  significant = "pval",
  p_cutoff = 0.05
)

cat("  ORA summary (by class):\n")
print(summary(ora_class_result))

png(file.path(base_dir, "03.Enrichment/ORA/ORA_by_class_barplot.png"),
  width = 1200, height = 1000, res = 150
)
print(ora_class_result$static_barPlot)
dev.off()

cat("  ✓ ORA by class complete\n\n")

# ORA by lipid category
cat("Running ORA by category...\n")
ora_category_result <- enrichment_ora(
  deSp_se,
  char = "Category",
  significant = "pval",
  p_cutoff = 0.05
)

cat("  ORA summary (by category):\n")
print(summary(ora_category_result))

png(file.path(base_dir, "03.Enrichment/ORA/ORA_by_category_barplot.png"),
  width = 1200, height = 1000, res = 150
)
print(ora_category_result$static_barPlot)
dev.off()

cat("  ✓ ORA by category complete\n")


##########################################################################################
# ANALYSIS COMPLETE
##########################################################################################

cat("\n========================================\n")
cat("POST-HOC ANALYSIS MILD vs SEVERE COMPLETED!\n")
cat("========================================\n")
