##########################################################################################
# PROFILING - POST-HOC COMPARISONS
##########################################################################################
#
# This script contains complete profiling analyses for both post-hoc comparisons:
# 1. HEALTHY vs MILD
# 2. MILD vs SEVERE
#
# REQUIREMENTS:
# - Variables must be defined before running this script:
#   * base_dir
#   * se_healthy_mild, processed_se_healthy_mild
#   * se_mild_severe, processed_se_mild_severe
#
##########################################################################################


##########################################################################################
# COMPARISON 1: HEALTHY vs MILD - COMPLETE PROFILING
##########################################################################################

cat("\n========================================\n")
cat("PROFILING: HEALTHY vs MILD\n")
cat("========================================\n")

# DATA PROCESSING QUALITY CHECKS
cat("\n--- Data Processing Quality ---\n")
data_process_plots_hm <- plot_data_process(se_healthy_mild, processed_se_healthy_mild)
summary(data_process_plots_hm)

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/00.Data_quality/BoxPlot_before_process.png"),
    width = 1400, height = 1000, res = 150)
print(data_process_plots_hm$static_boxPlot_before)
dev.off()

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/00.Data_quality/BoxPlot_after_process.png"),
    width = 1400, height = 1000, res = 150)
print(data_process_plots_hm$static_boxPlot_after)
dev.off()

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/00.Data_quality/DensityPlot_before_process.png"),
    width = 1400, height = 1000, res = 150)
print(data_process_plots_hm$static_densityPlot_before)
dev.off()

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/00.Data_quality/DensityPlot_after_process.png"),
    width = 1400, height = 1000, res = 150)
print(data_process_plots_hm$static_densityPlot_after)
dev.off()


# CROSS-SAMPLE VARIABILITY
cat("\n--- Cross-Sample Variability ---\n")
result_hm <- cross_sample_variability(se_healthy_mild)
summary(result_hm)

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/01.Cross-sample_variability/Expressed_lipid_numbers.png"),
    width = 1600, height = 1000, res = 150)
print(result_hm$static_lipid_number_barPlot)
dev.off()

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/01.Cross-sample_variability/Histogram_lipid_amount.png"),
    width = 1400, height = 1000, res = 150)
print(result_hm$static_lipid_amount_barPlot)
dev.off()

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/01.Cross-sample_variability/Lipid_abundance_distribution.png"),
    width = 1400, height = 1000, res = 150)
print(result_hm$static_lipid_distribution)
dev.off()


# DIMENSIONALITY REDUCTION
cat("\n--- Dimensionality Reduction ---\n")

## PCA
cat("Running PCA...\n")
result_pca_hm <- dr_pca(
  processed_se_healthy_mild, scaling = TRUE, centering = TRUE, clustering = "kmeans",
  cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
  eps = NULL, minPts = NULL, feature_contrib_pc = c(1, 2), plot_topN = 10)

summary(result_pca_hm)

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/02.Dimensionality_reduction/PCA/PCA_plot.png"),
    width = 1400, height = 1000, res = 150)
print(result_pca_hm$static_pca)
dev.off()

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/02.Dimensionality_reduction/PCA/Scree_plot.png"),
    width = 1400, height = 1000, res = 150)
print(result_pca_hm$static_screePlot)
dev.off()

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/02.Dimensionality_reduction/PCA/Feature_contribution.png"),
    width = 1400, height = 1000, res = 150)
print(result_pca_hm$static_feature_contribution)
dev.off()

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/02.Dimensionality_reduction/PCA/Variable_correlation.png"),
    width = 1400, height = 1000, res = 150)
print(result_pca_hm$static_variablePlot)
dev.off()

write_tsv(result_pca_hm$table_pca_contribution,
          file = file.path(base_dir, "01.Profiling/Healthy_vs_Mild/02.Dimensionality_reduction/PCA/table_pca_contribution.tsv"))
write_tsv(result_pca_hm$pca_rotated_data,
          file = file.path(base_dir, "01.Profiling/Healthy_vs_Mild/02.Dimensionality_reduction/PCA/table_pca_rotated_data.tsv"))


## t-SNE
cat("Running t-SNE...\n")
result_tsne_hm <- dr_tsne(
  processed_se_healthy_mild, pca = TRUE, perplexity = 5, max_iter = 500,
  clustering = "kmeans", cluster_num = 2, kmedoids_metric = NULL,
  distfun = NULL, hclustfun = NULL, eps = NULL, minPts = NULL)

summary(result_tsne_hm)

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/02.Dimensionality_reduction/t-SNE/t-SNE_plot.png"),
    width = 1400, height = 1000, res = 150)
print(result_tsne_hm$static_tsne)
dev.off()

write_tsv(result_tsne_hm$tsne_result,
          file = file.path(base_dir, "01.Profiling/Healthy_vs_Mild/02.Dimensionality_reduction/t-SNE/table_tsne_results.tsv"))


## UMAP
cat("Running UMAP...\n")
result_umap_hm <- dr_umap(
  processed_se_healthy_mild, n_neighbors = 15, scaling = TRUE,
  umap_metric = "euclidean", clustering = "kmeans", cluster_num = 2,
  kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
  eps = NULL, minPts = NULL)

summary(result_umap_hm)

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/02.Dimensionality_reduction/UMAP/UMAP_plot.png"),
    width = 1400, height = 1000, res = 150)
print(result_umap_hm$static_umap)
dev.off()

write_tsv(result_umap_hm$umap_result,
          file = file.path(base_dir, "01.Profiling/Healthy_vs_Mild/02.Dimensionality_reduction/UMAP/UMAP_results.tsv"))


# CORRELATION HEATMAPS
cat("\n--- Correlation Heatmaps ---\n")

## By sample
cat("Correlation by sample...\n")
result_heatmap_hm_sample <- heatmap_correlation(
  processed_se_healthy_mild, char = NULL, transform = "log10",
  correlation = "spearman", distfun = "euclidean",
  hclustfun = "ward.D2", type = "sample")

summary(result_heatmap_hm_sample)

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/03.Correlation_Heatmap/by_samples/Heatmap_by_samples.png"),
    width = 2000, height = 2000, res = 150)
print(result_heatmap_hm_sample$static_heatmap)
dev.off()


## By category
cat("Correlation by category...\n")
result_heatmap_hm_category <- heatmap_correlation(
  processed_se_healthy_mild, char = "Category", transform = "log10",
  correlation = "spearman", distfun = "euclidean",
  hclustfun = "ward.D2", type = "class")

summary(result_heatmap_hm_category)

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/03.Correlation_Heatmap/by_category/Heatmap_by_category.png"),
    width = 1500, height = 1500, res = 150)
print(result_heatmap_hm_category$static_heatmap)
dev.off()


## By class
cat("Correlation by lipid class...\n")
result_heatmap_hm_class <- heatmap_correlation(
  processed_se_healthy_mild, char = "class", transform = "log10",
  correlation = "spearman", distfun = "euclidean",
  hclustfun = "ward.D2", type = "class")

summary(result_heatmap_hm_class)

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/03.Correlation_Heatmap/by_class/Heatmap_by_class.png"),
    width = 1500, height = 1500, res = 150)
print(result_heatmap_hm_class$static_heatmap)
dev.off()

write_tsv(as.data.frame(result_heatmap_hm_class$corr_coef_matrix),
          file = file.path(base_dir, "01.Profiling/Healthy_vs_Mild/03.Correlation_Heatmap/by_class/Correlation_matrix.tsv"))


# LIPID CHARACTERISTICS PER SAMPLE
cat("\n--- Lipid Characteristics ---\n")

## By class
cat("Lipid profiling by class...\n")
result_lipid_hm <- lipid_profiling(processed_se_healthy_mild, char = "class")
summary(result_lipid_hm)

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/04.Lipid_characteristics/Lipid_class_per_sample.png"),
    width = 1000, height = 1000, res = 150)
print(result_lipid_hm$static_char_barPlot)
dev.off()

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/04.Lipid_characteristics/Lipid_composition_per_sample.png"),
    width = 1000, height = 1000, res = 150)
print(result_lipid_hm$static_lipid_composition)
dev.off()


## By Total.C
cat("Lipid profiling by Total.C...\n")
totalC_profile_hm <- lipid_profiling(processed_se_healthy_mild, char = "Total.C")

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/04.Lipid_characteristics/Total_C_per_sample.png"),
    width = 1000, height = 1000, res = 150)
print(totalC_profile_hm$static_char_barPlot)
dev.off()


## By Total.DB
cat("Lipid profiling by Total.DB...\n")
totalDB_profile_hm <- lipid_profiling(processed_se_healthy_mild, char = "Total.DB")

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/04.Lipid_characteristics/Total_DB_per_sample.png"),
    width = 1000, height = 1000, res = 150)
print(totalDB_profile_hm$static_char_barPlot)
dev.off()

cat("\n✓ Healthy vs Mild profiling completed!\n")


##########################################################################################
# COMPARISON 2: MILD vs SEVERE - COMPLETE PROFILING
##########################################################################################

cat("\n========================================\n")
cat("PROFILING: MILD vs SEVERE\n")
cat("========================================\n")

# DATA PROCESSING QUALITY CHECKS
cat("\n--- Data Processing Quality ---\n")
data_process_plots_ms <- plot_data_process(se_mild_severe, processed_se_mild_severe)
summary(data_process_plots_ms)

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/00.Data_quality/BoxPlot_before_process.png"),
    width = 1400, height = 1000, res = 150)
print(data_process_plots_ms$static_boxPlot_before)
dev.off()

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/00.Data_quality/BoxPlot_after_process.png"),
    width = 1400, height = 1000, res = 150)
print(data_process_plots_ms$static_boxPlot_after)
dev.off()

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/00.Data_quality/DensityPlot_before_process.png"),
    width = 1400, height = 1000, res = 150)
print(data_process_plots_ms$static_densityPlot_before)
dev.off()

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/00.Data_quality/DensityPlot_after_process.png"),
    width = 1400, height = 1000, res = 150)
print(data_process_plots_ms$static_densityPlot_after)
dev.off()


# CROSS-SAMPLE VARIABILITY
cat("\n--- Cross-Sample Variability ---\n")
result_ms <- cross_sample_variability(se_mild_severe)
summary(result_ms)

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/01.Cross-sample_variability/Expressed_lipid_numbers.png"),
    width = 1600, height = 1000, res = 150)
print(result_ms$static_lipid_number_barPlot)
dev.off()

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/01.Cross-sample_variability/Histogram_lipid_amount.png"),
    width = 1400, height = 1000, res = 150)
print(result_ms$static_lipid_amount_barPlot)
dev.off()

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/01.Cross-sample_variability/Lipid_abundance_distribution.png"),
    width = 1400, height = 1000, res = 150)
print(result_ms$static_lipid_distribution)
dev.off()


# DIMENSIONALITY REDUCTION
cat("\n--- Dimensionality Reduction ---\n")

## PCA
cat("Running PCA...\n")
result_pca_ms <- dr_pca(
  processed_se_mild_severe, scaling = TRUE, centering = TRUE, clustering = "kmeans",
  cluster_num = 2, kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
  eps = NULL, minPts = NULL, feature_contrib_pc = c(1, 2), plot_topN = 10)

summary(result_pca_ms)

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/02.Dimensionality_reduction/PCA/PCA_plot.png"),
    width = 1400, height = 1000, res = 150)
print(result_pca_ms$static_pca)
dev.off()

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/02.Dimensionality_reduction/PCA/Scree_plot.png"),
    width = 1400, height = 1000, res = 150)
print(result_pca_ms$static_screePlot)
dev.off()

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/02.Dimensionality_reduction/PCA/Feature_contribution.png"),
    width = 1400, height = 1000, res = 150)
print(result_pca_ms$static_feature_contribution)
dev.off()

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/02.Dimensionality_reduction/PCA/Variable_correlation.png"),
    width = 1400, height = 1000, res = 150)
print(result_pca_ms$static_variablePlot)
dev.off()

write_tsv(result_pca_ms$table_pca_contribution,
          file = file.path(base_dir, "01.Profiling/Mild_vs_Severe/02.Dimensionality_reduction/PCA/table_pca_contribution.tsv"))
write_tsv(result_pca_ms$pca_rotated_data,
          file = file.path(base_dir, "01.Profiling/Mild_vs_Severe/02.Dimensionality_reduction/PCA/table_pca_rotated_data.tsv"))


## t-SNE
cat("Running t-SNE...\n")
result_tsne_ms <- dr_tsne(
  processed_se_mild_severe, pca = TRUE, perplexity = 5, max_iter = 500,
  clustering = "kmeans", cluster_num = 2, kmedoids_metric = NULL,
  distfun = NULL, hclustfun = NULL, eps = NULL, minPts = NULL)

summary(result_tsne_ms)

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/02.Dimensionality_reduction/t-SNE/t-SNE_plot.png"),
    width = 1400, height = 1000, res = 150)
print(result_tsne_ms$static_tsne)
dev.off()

write_tsv(result_tsne_ms$tsne_result,
          file = file.path(base_dir, "01.Profiling/Mild_vs_Severe/02.Dimensionality_reduction/t-SNE/table_tsne_results.tsv"))


## UMAP
cat("Running UMAP...\n")
result_umap_ms <- dr_umap(
  processed_se_mild_severe, n_neighbors = 15, scaling = TRUE,
  umap_metric = "euclidean", clustering = "kmeans", cluster_num = 2,
  kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL,
  eps = NULL, minPts = NULL)

summary(result_umap_ms)

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/02.Dimensionality_reduction/UMAP/UMAP_plot.png"),
    width = 1400, height = 1000, res = 150)
print(result_umap_ms$static_umap)
dev.off()

write_tsv(result_umap_ms$umap_result,
          file = file.path(base_dir, "01.Profiling/Mild_vs_Severe/02.Dimensionality_reduction/UMAP/UMAP_results.tsv"))


# CORRELATION HEATMAPS
cat("\n--- Correlation Heatmaps ---\n")

## By sample
cat("Correlation by sample...\n")
result_heatmap_ms_sample <- heatmap_correlation(
  processed_se_mild_severe, char = NULL, transform = "log10",
  correlation = "spearman", distfun = "euclidean",
  hclustfun = "ward.D2", type = "sample")

summary(result_heatmap_ms_sample)

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/03.Correlation_Heatmap/by_samples/Heatmap_by_samples.png"),
    width = 2000, height = 2000, res = 150)
print(result_heatmap_ms_sample$static_heatmap)
dev.off()


## By category
cat("Correlation by category...\n")
result_heatmap_ms_category <- heatmap_correlation(
  processed_se_mild_severe, char = "Category", transform = "log10",
  correlation = "spearman", distfun = "euclidean",
  hclustfun = "ward.D2", type = "class")

summary(result_heatmap_ms_category)

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/03.Correlation_Heatmap/by_category/Heatmap_by_category.png"),
    width = 1500, height = 1500, res = 150)
print(result_heatmap_ms_category$static_heatmap)
dev.off()


## By class
cat("Correlation by lipid class...\n")
result_heatmap_ms_class <- heatmap_correlation(
  processed_se_mild_severe, char = "class", transform = "log10",
  correlation = "spearman", distfun = "euclidean",
  hclustfun = "ward.D2", type = "class")

summary(result_heatmap_ms_class)

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/03.Correlation_Heatmap/by_class/Heatmap_by_class.png"),
    width = 1500, height = 1500, res = 150)
print(result_heatmap_ms_class$static_heatmap)
dev.off()

write_tsv(as.data.frame(result_heatmap_ms_class$corr_coef_matrix),
          file = file.path(base_dir, "01.Profiling/Mild_vs_Severe/03.Correlation_Heatmap/by_class/Correlation_matrix.tsv"))


# LIPID CHARACTERISTICS PER SAMPLE
cat("\n--- Lipid Characteristics ---\n")

## By class
cat("Lipid profiling by class...\n")
result_lipid_ms <- lipid_profiling(processed_se_mild_severe, char = "class")
summary(result_lipid_ms)

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/04.Lipid_characteristics/Lipid_class_per_sample.png"),
    width = 1000, height = 1000, res = 150)
print(result_lipid_ms$static_char_barPlot)
dev.off()

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/04.Lipid_characteristics/Lipid_composition_per_sample.png"),
    width = 1000, height = 1000, res = 150)
print(result_lipid_ms$static_lipid_composition)
dev.off()


## By Total.C
cat("Lipid profiling by Total.C...\n")
totalC_profile_ms <- lipid_profiling(processed_se_mild_severe, char = "Total.C")

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/04.Lipid_characteristics/Total_C_per_sample.png"),
    width = 1000, height = 1000, res = 150)
print(totalC_profile_ms$static_char_barPlot)
dev.off()


## By Total.DB
cat("Lipid profiling by Total.DB...\n")
totalDB_profile_ms <- lipid_profiling(processed_se_mild_severe, char = "Total.DB")

png(file.path(base_dir, "01.Profiling/Mild_vs_Severe/04.Lipid_characteristics/Total_DB_per_sample.png"),
    width = 1000, height = 1000, res = 150)
print(totalDB_profile_ms$static_char_barPlot)
dev.off()

cat("\n✓ Mild vs Severe profiling completed!\n")

cat("\n========================================\n")
cat("ALL POST-HOC PROFILING COMPLETED!\n")
cat("========================================\n")

