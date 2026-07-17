##########################################################################################
# LIPIDOMICS ANALYSIS: HEALTHY vs MILD vs SEVERE PATIENTS
##########################################################################################
#
# ANALYSIS STRATEGY:
# ------------------
# This script performs POST-HOC differential expression analysis using pairwise t-tests
# instead of One-way ANOVA. This approach is MORE POWERFUL for detecting differences
# between specific groups.
#
# COMPARISONS PERFORMED:
# 1. HEALTHY vs MILD: Identifies lipids altered in mild disease compared to healthy controls
# 2. MILD vs SEVERE: Identifies lipids that change with disease severity
#
# ADVANTAGES OF POST-HOC T-TESTS OVER ANOVA:
# - Higher statistical power for pairwise comparisons
# - More sensitive detection of specific group differences
# - Better suited for unbalanced group sizes
# - Provides directional information (up/down regulation)
# - Enables separate volcano plots and MA plots for each comparison
#
# DATA ORGANIZATION:
# - Three filtered datasets are created: full (all 3 groups), healthy+mild, mild+severe
# - All profiling and PCA uses the full dataset (3 groups)
# - Differential expression uses the filtered pairwise datasets
#
##########################################################################################

##########################
# SETUP
##########################

# Step 1: Install devtools
install.packages("devtools")

# Step 2: Install BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Step 3: Install LipidSigR
## Update repositories
options(repos = c(
  CRAN = "https://cloud.r-project.org/",
  BiocManager::repositories()))

## Install dependencies and package
devtools::install_github(
  "BioinfOMICS/LipidSigR",
  build_vignettes = TRUE, dependencies = TRUE)

devtools::install_github("ricoderks/Rcpm")

# Install dependencies:
BiocManager::install(
  c('impute', 'fgsea', 'gatom', 'mixOmics', 'S4Vectors', 'BiocGenerics',
    'SummarizedExperiment', 'rgoslin'))

install.packages(
  c('', 'magrittr', 'plotly', 'tidyverse', 'factoextra', 'ggthemes',
    'ggforce', 'Hmisc', 'hwordcloud', 'heatmaply', 'iheatmapr', 'Rtsne', 'uwot',
    'wordcloud', 'rsample', 'ranger', 'caret', 'yardstick', 'fastshap',
    'SHAPforxgboost', 'visNetwork', 'tidygraph', 'ggraph'))

devtools::install_github("ctlab/mwcsr")

BiocManager::install("rgoslin")

##########################
# LOAD DATA
##########################

# Read files
library(readr)
library(LipidSigR)
#ls("package:LipidSigR")

group_info_threeGroup <- read_tsv("/Users/loictalignani/research/project/lipidomics/data/lipidsig_datasets/healthy_vs_sick_patients/group_information_table_healthy_vs_sick_patients_D0.tsv")
head(group_info_threeGroup, 5)

abundance_threeGroup <- read_tsv("/Users/loictalignani/research/project/lipidomics/data/lipidsig_datasets/healthy_vs_sick_patients/healthy_sick_lipidomics.tsv")
head(abundance_threeGroup)

# map lipid characteristics by rgoslin
library(rgoslin)
library(dplyr)

parse_lipid <- rgoslin::parseLipidNames(lipidNames=abundance_threeGroup$feature)

# filter lipid recognized by rgoslin
recognized_lipid <- parse_lipid$Original.Name[
  which(parse_lipid$Grammar != 'NOT_PARSEABLE')]
abundance <- abundance_threeGroup %>% 
  dplyr::filter(feature %in% recognized_lipid)
goslin_annotation <- parse_lipid %>% 
  dplyr::filter(Original.Name %in% recognized_lipid)

head(abundance[, 1:6], 5)

head(goslin_annotation[, 1:6], 5)

# Construct SE object for all three groups (for profiling and PCA)
se <- as_summarized_experiment(
  abundance, goslin_annotation, group_info=group_info_threeGroup,
  se_type='de_multiple', paired_sample=NULL)

se@assays@data@listData

# data processing
processed_se <- data_process(
  se, exclude_missing=TRUE, exclude_missing_pct=70,
  replace_na_method='min', replace_na_method_ref=0.5,
  normalization='Percentage', transform='log10')

processed_se@metadata$processed_abund

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

####################################################
# PREPARE FILTERED DATASETS FOR POST-HOC COMPARISONS
######################################################

##########################
# COMPARISON 1: HEALTHY vs MILD
##########################

# Load group information for healthy vs mild
group_info_healthy_mild <- read_tsv("/Users/loictalignani/research/project/lipidomics/data/lipidsig_datasets/healthy_vs_sick_patients/group_information_table_healthy_vs_MILD_D0.tsv")
head(group_info_healthy_mild, 5)

# Load abundance data for healthy vs mild
abundance_healthy_mild <- read_tsv("/Users/loictalignani/research/project/lipidomics/data/lipidsig_datasets/healthy_vs_sick_patients/healthy_MILD_lipidomics.tsv")
head(abundance_healthy_mild)

# Parse lipid names with rgoslin for healthy vs mild
parse_lipid_hm <- rgoslin::parseLipidNames(lipidNames=abundance_healthy_mild$feature)

# Filter lipid recognized by rgoslin
recognized_lipid_hm <- parse_lipid_hm$Original.Name[
  which(parse_lipid_hm$Grammar != 'NOT_PARSEABLE')]
abundance_hm <- abundance_healthy_mild %>%
  dplyr::filter(feature %in% recognized_lipid_hm)
goslin_annotation_hm <- parse_lipid_hm %>%
  dplyr::filter(Original.Name %in% recognized_lipid_hm)

# Construct SE object for healthy vs mild
se_healthy_mild <- as_summarized_experiment(
  abundance_hm, goslin_annotation_hm, group_info=group_info_healthy_mild,
  se_type='de_two', paired_sample=FALSE)

# Data processing for healthy vs mild
processed_se_healthy_mild <- data_process(
  se_healthy_mild, exclude_missing=TRUE, exclude_missing_pct=50,
  replace_na_method='min', replace_na_method_ref=0.5,
  normalization='Percentage', transform='log10')

##########################
# COMPARISON 2: MILD vs SEVERE
##########################

# Load group information for mild vs severe
group_info_mild_severe <- read_tsv("/Users/loictalignani/research/project/lipidomics/data/lipidsig_datasets/healthy_vs_sick_patients/group_information_table_MILD_vs_SEVERE_patients_D0.tsv")
head(group_info_mild_severe, 5)

# Load abundance data for mild vs severe
abundance_mild_severe <- read_tsv("/Users/loictalignani/research/project/lipidomics/data/lipidsig_datasets/healthy_vs_sick_patients/MILD_SEVERE_lipidomics.tsv")
head(abundance_mild_severe)

# Parse lipid names with rgoslin for mild vs severe
parse_lipid_ms <- rgoslin::parseLipidNames(lipidNames=abundance_mild_severe$feature)

# Filter lipid recognized by rgoslin
recognized_lipid_ms <- parse_lipid_ms$Original.Name[
  which(parse_lipid_ms$Grammar != 'NOT_PARSEABLE')]
abundance_ms <- abundance_mild_severe %>%
  dplyr::filter(feature %in% recognized_lipid_ms)
goslin_annotation_ms <- parse_lipid_ms %>%
  dplyr::filter(Original.Name %in% recognized_lipid_ms)

# Construct SE object for mild vs severe
se_mild_severe <- as_summarized_experiment(
  abundance_ms, goslin_annotation_ms, group_info=group_info_mild_severe,
  se_type='de_two', paired_sample=FALSE)

# Data processing for mild vs severe
processed_se_mild_severe <- data_process(
  se_mild_severe, exclude_missing=TRUE, exclude_missing_pct=50,
  replace_na_method='min', replace_na_method_ref=0.5,
  normalization='Percentage', transform='log10')


##########################################################################################
# PROFILING - POST-HOC COMPARISONS
##########################################################################################
#
# This section performs profiling for both post-hoc comparisons:
# 1. HEALTHY vs MILD
# 2. MILD vs SEVERE
#
# For each comparison, we perform:
# - Data processing quality checks (before/after plots)
# - Cross-sample variability analysis
# - Dimensionality reduction (PCA, t-SNE, UMAP)
# - Correlation heatmaps
# - Lipid characteristics per sample
##########################################################################################

# CREATE DIRECTORY STRUCTURE
cat("\n=== CREATING DIRECTORY STRUCTURE ===\n")

# Base directories for both comparisons
base_dir <- "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0"

# Healthy vs Mild directories
dir.create(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/00.Data_quality"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/01.Cross-sample_variability"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/02.Dimensionality_reduction/PCA"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/02.Dimensionality_reduction/t-SNE"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/02.Dimensionality_reduction/UMAP"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/03.Correlation_Heatmap/by_samples"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/03.Correlation_Heatmap/by_category"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/03.Correlation_Heatmap/by_class"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/04.Lipid_characteristics"), recursive=TRUE, showWarnings=FALSE)

# Mild vs Severe directories
dir.create(file.path(base_dir, "01.Profiling/Mild_vs_Severe/00.Data_quality"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(base_dir, "01.Profiling/Mild_vs_Severe/01.Cross-sample_variability"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(base_dir, "01.Profiling/Mild_vs_Severe/02.Dimensionality_reduction/PCA"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(base_dir, "01.Profiling/Mild_vs_Severe/02.Dimensionality_reduction/t-SNE"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(base_dir, "01.Profiling/Mild_vs_Severe/02.Dimensionality_reduction/UMAP"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(base_dir, "01.Profiling/Mild_vs_Severe/03.Correlation_Heatmap/by_samples"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(base_dir, "01.Profiling/Mild_vs_Severe/03.Correlation_Heatmap/by_category"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(base_dir, "01.Profiling/Mild_vs_Severe/03.Correlation_Heatmap/by_class"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(base_dir, "01.Profiling/Mild_vs_Severe/04.Lipid_characteristics"), recursive=TRUE, showWarnings=FALSE)

cat("Directory structure created successfully!\n")


##########################################################################################
# COMPARISON 1: HEALTHY vs MILD - PROFILING
##########################################################################################

cat("\n=== PROFILING: HEALTHY vs MILD ===\n")

# DATA PROCESSING QUALITY CHECKS
cat("\n--- Data Processing Quality ---\n")
data_process_plots_hm <- plot_data_process(se_healthy_mild, processed_se_healthy_mild)

summary(data_process_plots_hm)

# Save plots: Box plot before/after data processing
png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/00.Data_quality/BoxPlot_before_process.png"),
    width = 1400, height = 1000, res = 150)
print(data_process_plots_hm$static_boxPlot_before)
dev.off()

png(file.path(base_dir, "01.Profiling/Healthy_vs_Mild/00.Data_quality/BoxPlot_after_process.png"),
    width = 1400, height = 1000, res = 150)
print(data_process_plots_hm$static_boxPlot_after)
dev.off()

# Save plots: Density plot before/after data processing
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

# Save plots
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
  processed_se_healthy_mild, scaling=TRUE, centering=TRUE, clustering='kmeans',
  cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL,
  eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10)

# result summary
summary(result_pca)

# view result: PCA plot
result_pca$interactive_pca
static_pca <- result_pca$static_pca
png(
  "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/01.Profiling/02.Results_D0/01.Dimensionality_reduction/PCA/static_PCA_plot.png",
  width = 1400, height = 1000, res = 150
)
print(static_pca)
dev.off()

# view result: scree plot of top 10 principle components
result_pca$static_screePlot
static_screePlot <- result_pca$static_screePlot
png(
  "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/01.Profiling/02.Results_D0/01.Dimensionality_reduction/PCA/Static_explained_variance.png",
  width = 1400, height = 1000, res = 150
)
print(static_screePlot)
dev.off()

# view result: correlation circle plot of PCA variables
result_pca$static_feature_contribution
static_feature_contribution <- result_pca$static_feature_contribution
png(
  "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/01.Profiling/02.Results_D0/01.Dimensionality_reduction/PCA/static_feature_contribution.png",
  width = 1400, height = 1000, res = 150
)
print(static_feature_contribution)
dev.off()

# view result: Correlation of contribution of top 10 features
result_pca$static_variablePlot
static_variablePlot <- result_pca$static_variablePlot
png(
  "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/01.Profiling/02.Results_D0/01.Dimensionality_reduction/PCA/corr_circle_plot.png",
  width = 1400, height = 1000, res = 150
)
print(static_variablePlot)
dev.off()

# save results
write_tsv(result_pca$table_pca_contribution, file = "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/01.Profiling/02.Results_D0/01.Dimensionality_reduction/PCA/table_pca_contribution.tsv")
write_tsv(result_pca$pca_rotated_data, file = "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/01.Profiling/02.Results_D0/01.Dimensionality_reduction/PCA/table_pca_rotated_data.tsv")
# T-SNE


# conduct t-SNE
result_tsne <- dr_tsne(
  processed_se, pca=TRUE, perplexity=5, max_iter=500, clustering='kmeans',
  cluster_num=3, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, 
  eps=NULL, minPts=NULL)

# result summary
summary(result_tsne)

# view result: t-SNE plot
result_tsne$static_tsne
static_tsne <- result_tsne$static_tsne
png(
  "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/01.Profiling/02.Results_D0/01.Dimensionality_reduction/t-SNE/static_t-SNE_plot.png",
  width = 1400, height = 1000, res = 150
)
print(static_tsne)
dev.off()

# save results
write_tsv(result_tsne$tsne_result, file = "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/01.Profiling/02.Results_D0/01.Dimensionality_reduction/t-SNE/table_tsne_results.tsv")


# UMAP
# conduct UMAP
result_umap <- dr_umap(
  processed_se, n_neighbors=15, scaling=TRUE, umap_metric='euclidean',
  clustering='kmeans', cluster_num=3, kmedoids_metric=NULL,
  distfun=NULL, hclustfun=NULL, eps=NULL, minPts=NULL)

# result summary
summary(result_umap)

# view result: UMAP plot
result_umap$static_umap
static_umap <- result_umap$static_umap
png(
  "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/01.Profiling/02.Results_D0/01.Dimensionality_reduction/UMAP/static_UMAP_plot.png",
  width = 1400, height = 1000, res = 150
)
print(static_umap)
dev.off()

# save umap results
write_tsv(result_umap$umap_result, file = "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/01.Profiling/02.Results_D0/01.Dimensionality_reduction/UMAP/UMAP_results.tsv")

result_umap$interactive_umap




# CORRELATION HEATMAPS

# List all characteristics (char=...) available
list_lipid_char(processed_se)

# correlation calculation by sample
result_heatmap <- heatmap_correlation(
  processed_se, char=NULL, transform='log10', correlation='spearman', 
  distfun='euclidean', hclustfun='ward.D2', type='sample')

# result summary          
summary(result_heatmap)

result_heatmap$static_heatmap
static_heatmap <- result_heatmap$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/01.Profiling/02.Results_D0/02.Correlation_Heatmap/by_samples/Correlation_Heatmap_by_samples_spearman_Eucl_WD2.png",
  width = 2000, height = 2000, res = 150
)
print(static_heatmap)
dev.off()

# correlation calculation by category
result_heatmap <- heatmap_correlation(
  processed_se, char="Category", transform='log10', correlation='spearman', 
  distfun='euclidean', hclustfun='ward.D2', type='class')

# result summary          
summary(result_heatmap)

result_heatmap$static_heatmap
static_heatmap <- result_heatmap$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/01.Profiling/02.Results_D0/02.Correlation_Heatmap/by_category/Correlation_Heatmap_by_lipid_category_spearman_Eucl_WD2.png",
    width = 1500, height = 1500, res = 150
)
print(static_heatmap)
dev.off()

# correlation calculation by class and by sample
result_heatmap <- heatmap_correlation(
  processed_se, char="class", transform='log10', correlation='spearman', 
  distfun='euclidean', hclustfun='ward.D2', type='sample')

result_heatmap$static_heatmap


# correlation calculation by lipid class
result_heatmap <- heatmap_correlation(
  processed_se, char="class", transform='log10', correlation='spearman', 
  distfun='euclidean', hclustfun='ward.D2', type='class')

result_heatmap$static_heatmap

staticheatmap <- result_heatmap$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/01.Profiling/02.Results_D0/02.Correlation_Heatmap/by_class/Correlation_Heatmap_by_lipid_class_spearman_Eucl_WD2.png",
    width = 1500, height = 1500, res = 150
)
print(staticheatmap)
dev.off()
# Export de la matrice de corrélation en TSV
mat_corr <- result_heatmap$corr_coef_matrix
readr::write_tsv(as.data.frame(mat_corr),
                 file = "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/01.Profiling/02.Results_D0/02.Correlation_Heatmap/by_class/Correlation_Heatmap_by_lipid_class_spearman_Eucl_WD2.tsv")



# LIPID CHARACTERISTICS PER SAMPLE
# calculate lipid expression of selected characteristic
result_lipid <- lipid_profiling(processed_se, char="class")

# result summary
summary(result_lipid)

# view result: bar plot of class per sample
result_lipid$interactive_char_barPlot
static_char_barPlot <- result_lipid$static_char_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_expression/Lipid_class_per_sample.png",
    width = 1000, height = 1000, res = 150
)
print(static_char_barPlot)
dev.off()

# view result: stacked horizontal bar chart
result_lipid$interactive_lipid_composition
static_lipid_composition <- result_lipid$static_lipid_composition
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_expression/Lipid_composition_per_sample.png",
    width = 1000, height = 1000, res = 150
)
print(static_lipid_composition)
dev.off()


# calculate lipid expression of selected characteristic: Total.C
totalC_profile <- lipid_profiling(processed_se, char = "Total.C")
totalC_profile$interactive_char_barPlot
static_char_barPlot <- totalC_profile$static_char_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_expression/Total.C_per_sample.png",
    width = 1000, height = 1000, res = 150
)
print(static_char_barPlot)
dev.off()


# calculate lipid expression of selected characteristic: Total.DB
totalDB_profile <- lipid_profiling(processed_se, char = "Total.DB")
totalDB_profile$interactive_char_barPlot
static_char_barPlot <- totalDB_profile$static_char_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_expression/Total.DB_per_sample.png",
    width = 1000, height = 1000, res = 150
)
print(static_char_barPlot)
dev.off()

# calculate lipid expression of selected characteristic: Total.OH
totalOH_profile <- lipid_profiling(processed_se, char = "Total.OH")
totalOH_profile$interactive_char_barPlot
static_char_barPlot <- totalOH_profile$static_char_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_expression/Total.OH_per_sample.png",
    width = 1000, height = 1000, res = 150
)
print(static_char_barPlot)
dev.off()



# calculate lipid expression of selected characteristic: Total.FA
totalFA_profile <- lipid_profiling(processed_se, char = "Total.FA")
totalFA_profile$interactive_char_barPlot
static_char_barPlot <- totalFA_profile$static_char_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_expression/Total.FA_per_sample.png",
    width = 3000, height = 3000, res = 150
)
print(static_char_barPlot)
dev.off()

##################################################
# POST-HOC DIFFERENTIAL EXPRESSION ANALYSIS
# Two-by-two comparisons with t-test (more powerful than ANOVA)
##################################################

################################
# COMPARISON 1: HEALTHY vs MILD
################################

# STEP 1: CONDUCT DIFFERENTIAL EXPRESSION
cat("\n=== DIFFERENTIAL EXPRESSION ANALYSIS: HEALTHY vs MILD ===\n")

deSp_se_healthy_mild <- deSp_twoGroup(
  processed_se_healthy_mild, ref_group='healthy', test='t-test',
  significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')

# extract results in SE
res_list_healthy_mild <- extract_summarized_experiment(deSp_se_healthy_mild)
# summary of extract results
cat("\nSummary of results (healthy vs mild):\n")
summary(res_list_healthy_mild)

# plot differential expression analysis result
deSp_plot_healthy_mild <- plot_deSp_twoGroup(deSp_se_healthy_mild)

# result summary
summary(deSp_plot_healthy_mild)

# view result: lollipop chart
deSp_plot_healthy_mild$static_de_lipid
static_de_lipid_hm <- deSp_plot_healthy_mild$static_de_lipid
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Healthy_vs_Mild/Lollipop_chart/lollipop_chart.png",
    width = 1000, height = 1000, res = 150
)
print(static_de_lipid_hm)
dev.off()

# view result: MA plot
deSp_plot_healthy_mild$interactive_maPlot
static_maPlot_hm <- deSp_plot_healthy_mild$static_maPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Healthy_vs_Mild/MA_plot/MA_plot.png",
    width = 1000, height = 1000, res = 150
)
print(static_maPlot_hm)
dev.off()

# view result: Volcano plot
deSp_plot_healthy_mild$interactive_volcanoPlot
static_volcanoPlot_hm <- deSp_plot_healthy_mild$static_volcanoPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Healthy_vs_Mild/Volcano_plot/Volcano_plot.png",
    width = 1000, height = 1000, res = 150
)
print(static_volcanoPlot_hm)
dev.off()

################################
# COMPARISON 2: MILD vs SEVERE
################################

# STEP 1: CONDUCT DIFFERENTIAL EXPRESSION
cat("\n=== DIFFERENTIAL EXPRESSION ANALYSIS: MILD vs SEVERE ===\n")

deSp_se_mild_severe <- deSp_twoGroup(
  processed_se_mild_severe, ref_group='mild', test='t-test',
  significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')

# extract results in SE
res_list_mild_severe <- extract_summarized_experiment(deSp_se_mild_severe)
# summary of extract results
cat("\nSummary of results (mild vs severe):\n")
summary(res_list_mild_severe)

# plot differential expression analysis result
deSp_plot_mild_severe <- plot_deSp_twoGroup(deSp_se_mild_severe)

# result summary
summary(deSp_plot_mild_severe)

# view result: lollipop chart
deSp_plot_mild_severe$static_de_lipid
static_de_lipid_ms <- deSp_plot_mild_severe$static_de_lipid
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Mild_vs_Severe/Lollipop_chart/lollipop_chart.png",
    width = 1000, height = 1000, res = 150
)
print(static_de_lipid_ms)
dev.off()

# view result: MA plot
deSp_plot_mild_severe$interactive_maPlot
static_maPlot_ms <- deSp_plot_mild_severe$static_maPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Mild_vs_Severe/MA_plot/MA_plot.png",
    width = 1000, height = 1000, res = 150
)
print(static_maPlot_ms)
dev.off()

# view result: Volcano plot
deSp_plot_mild_severe$interactive_volcanoPlot
static_volcanoPlot_ms <- deSp_plot_mild_severe$static_volcanoPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Mild_vs_Severe/Volcano_plot/Volcano_plot.png",
    width = 1000, height = 1000, res = 150
)
print(static_volcanoPlot_ms)
dev.off()

# Export results tables
write_tsv(deSp_plot_healthy_mild$table_de_lipid, file = "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Healthy_vs_Mild/table_diffexp.tsv")
write_tsv(deSp_plot_mild_severe$table_de_lipid, file = "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Mild_vs_Severe/table_diffexp.tsv")

##################################################
# INDIVIDUAL LIPID ABUNDANCE BOXPLOTS
# For specific lipids of interest from differential expression
##################################################

# Example: plot abundance box plot for a lipid of interest
# You can identify lipids of interest from the differential expression results above
# and create individual boxplots using boxPlot_feature_twoGroup()

# Example for healthy vs mild comparison:
# boxPlot_result_hm <- boxPlot_feature_twoGroup(
#   processed_se_healthy_mild, feature='LPC 16:0', ref_group='healthy',
#   test='t-test', transform='log10')
# boxPlot_result_hm$static_boxPlot
# png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Healthy_vs_Mild/Boxplots/LPC_16_0.png",
#     width = 1000, height = 1000, res = 150)
# print(boxPlot_result_hm$static_boxPlot)
# dev.off()

# Example for mild vs severe comparison:
# boxPlot_result_ms <- boxPlot_feature_twoGroup(
#   processed_se_mild_severe, feature='LPC 16:0', ref_group='mild',
#   test='t-test', transform='log10')
# boxPlot_result_ms$static_boxPlot
# png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Mild_vs_Severe/Boxplots/LPC_16_0.png",
#     width = 1000, height = 1000, res = 150)
# print(boxPlot_result_ms$static_boxPlot)
# dev.off()

##################################################
# LIPID CHARACTERISTICS DIFFERENTIAL EXPRESSION ANALYSIS
# Post-hoc comparisons for lipid characteristics
##################################################

# Check available characteristics
list_lipid_char(processed_se)

################################
# COMPARISON 1: HEALTHY vs MILD - Characteristics
################################

# conduct differential expression of lipid characteristics for "Total.C"
deChar_se_Total_C_hm <- deChar_twoGroup(
  processed_se_healthy_mild, char="Total.C", ref_group="healthy", test='t-test', 
  post_hoc_sig = "padj", post_hoc_p_cutoff = 0.05, transform='log10')

# plot differential expression analysis results
deChar_plot <- plot_deChar_multiGroup(deChar_se)

# result summary
summary(deChar_plot)

# view result: bar plot of selected `Total.C`
deChar_plot$interactive_barPlot
static_barPlot <- deChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/TotalC/Barplot_of_characteristics_analysis_TotalC.png",
    width = 2000, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# conduct differential expression of lipid characteristics for Fatty acid properties "Total.DB"
deChar_se <- deChar_multiGroup(
  processed_se, char="Total.DB", ref_group="healthy", post_hoc ='One-way ANOVA', 
  post_hoc_sig = "padj", post_hoc_p_cutoff = 0.05, transform='log10')

# plot differential expression analysis results
deChar_plot <- plot_deChar_multiGroup(deChar_se)

# result summary
summary(deChar_plot)

# view result: bar plot of selected `Total.DB`
deChar_plot$interactive_barPlot
static_barPlot <- deChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/TotalDB/Barplot_of_characteristics_analysis_TotalDB.png",
    width = 2000, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# conduct differential expression of lipid characteristics for Fatty acid properties "Total.FA"
deChar_se <- deChar_multiGroup(
  processed_se, char="Total.FA", ref_group="healthy", post_hoc ='One-way ANOVA', 
  post_hoc_sig = "padj", post_hoc_p_cutoff = 0.05, transform='log10')

# plot differential expression analysis results
deChar_plot <- plot_deChar_multiGroup(deChar_se)

# result summary
summary(deChar_plot)

# view result: bar plot of selected `Total.FA`
deChar_plot$interactive_barPlot
static_barPlot <- deChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/TotalFA/Barplot_of_characteristics_analysis_TotalFA.png",
    width = 2000, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# conduct differential expression of lipid characteristics for Fatty acid properties "Total.OH"
deChar_se <- deChar_multiGroup(
  processed_se, char="Total.OH", ref_group="healthy", post_hoc ='One-way ANOVA', 
  post_hoc_sig = "padj", post_hoc_p_cutoff = 0.05, transform='log10')

# plot differential expression analysis results
deChar_plot <- plot_deChar_multiGroup(deChar_se)

# result summary
summary(deChar_plot)

# view result: bar plot of selected `Total.OH`
deChar_plot$interactive_barPlot
static_barPlot <- deChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/TotalOH/Barplot_of_characteristics_analysis_TotalOH.png",
    width = 2000, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# conduct differential expression of lipid characteristics for "class"
deChar_se <- deChar_multiGroup(
  processed_se, char="class", ref_group="healthy", post_hoc ='One-way ANOVA', 
  post_hoc_sig = "padj", post_hoc_p_cutoff = 0.05, transform='log10')

# plot differential expression analysis results
deChar_plot <- plot_deChar_multiGroup(deChar_se)

# result summary
summary(deChar_plot)

# view result: bar plot of selected `char`
deChar_plot$interactive_barPlot
static_barPlot <- deChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Class/Barplot_of_characteristics_analysis_class.png",
    width = 2000, height = 800, res = 150
    )
print(static_barPlot)
dev.off()

# conduct differential expression of lipid characteristics for "Category"
deChar_se <- deChar_multiGroup(
  processed_se, char="Category", ref_group="healthy", post_hoc ='One-way ANOVA', 
  post_hoc_sig = "padj", post_hoc_p_cutoff = 0.05, transform='log10')

# plot differential expression analysis results
deChar_plot <- plot_deChar_multiGroup(deChar_se)

# result summary
summary(deChar_plot)

# view result: bar plot of selected `char`
deChar_plot$interactive_barPlot
static_barPlot <- deChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Category/Barplot_of_characteristics_analysis_category.png",
    width = 2000, height = 800, res = 150
)
print(static_barPlot)
dev.off()

###################################
# SPECIES BASED AGGREGATION

#########
# TOTAL.C

# subgroup differential expression of lipid characteristics
subChar_se <- subChar_twoGroup(
  processed_se, char="Total.C", subChar="class", ref_group="mild", 
  test='t-test', significant="pval", p_cutoff=0.05, 
  FC_cutoff=1, transform='log10')

# get subChar_feature list
unique(
  extract_summarized_experiment(subChar_se)$all_deChar_result$sub_feature)

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="LPC O-")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_LPCO-.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="TG")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_TG.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PC")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_PC.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="CE")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_CE.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="DG")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_DG.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="LPC")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_LPC.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="FA")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_FA.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="CAR")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_CAR.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PI")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_PI.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="Cer")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_Cer.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PC O-")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_PC O-.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

###########
# TOTAL.DB

# subgroup differential expression of lipid characteristics
subChar_se <- subChar_twoGroup(
  processed_se, char="Total.DB", subChar="class", ref_group="mild", 
  test='t-test', significant="pval", p_cutoff=0.05, 
  FC_cutoff=1, transform='log10')

# get subChar_feature list
unique(
  extract_summarized_experiment(subChar_se)$all_deChar_result$sub_feature)

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="LPC O-")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_LPCO-.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="TG")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_TG.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PC")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_PC.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="CE")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_CE.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="DG")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_DG.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="LPC")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_LPC.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="FA")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_FA.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="CAR")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_CAR.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PI")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_PI.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="Cer")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_Cer.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PC O-")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_PC O-.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()


###########
# TOTAL.OH

# subgroup differential expression of lipid characteristics
subChar_se <- subChar_twoGroup(
  processed_se, char="Total.OH", subChar="class", ref_group="mild", 
  test='t-test', significant="pval", p_cutoff=0.05, 
  FC_cutoff=1, transform='log10')

# get subChar_feature list
unique(
  extract_summarized_experiment(subChar_se)$all_deChar_result$sub_feature)

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="LPC O-")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_LPCO-.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="TG")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_TG.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PC")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_PC.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="CE")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_CE.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="DG")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_DG.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="LPC")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_LPC.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="FA")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_FA.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="CAR")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_CAR.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PI")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_PI.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="Cer")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_Cer.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PC O-")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_PC O-.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()


###################################
# FATTY ACID BASED AGGREGATION

###############
# FA PER CATEGORY

list_lipid_char(processed_se)


# subgroup differential expression of lipid characteristics
subChar_se <- subChar_twoGroup(
  processed_se, char="FA.C", subChar="Category", ref_group="mild", 
  test='t-test', significant="pval", p_cutoff=0.05, 
  FC_cutoff=1, transform='log10')

# get subChar_feature list
unique(
  extract_summarized_experiment(subChar_se)$all_deChar_result$sub_feature)

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="Glycerophospholipids [GP]")

summary(subChar_plot)

# view result: bar plot of `subChar_feature`
subChar_plot$interactive_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Fatty_acids_based_aggregation/per_category/Glycerophospholipids.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="Sterol lipids [ST]")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Fatty_acids_based_aggregation/per_category/Sterol_lipids.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="Fatty acyls [FA]")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Fatty_acids_based_aggregation/per_category/Fatty_acyls.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()



#####################################
# STEP 2: DIMENSION REDUCTION: PLS-DA
# conduct PLSDA
result_plsda <- dr_plsda(
  deSp_se, ncomp=2, scaling=TRUE, clustering='group_info', cluster_num=3, 
  kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, eps=NULL, minPts=NULL)

# result summary
summary(result_plsda)
write.csv(result_plsda$plsda_result, file = "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/02.Dimensionality_reduction/PLS-DA/Table_of_sample_variate_D0.csv", row.names = FALSE)
write.csv(result_plsda$table_plsda_loading, file = "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/02.Dimensionality_reduction/PLS-DA/Table_of_sample_loading_D0.csv", row.names = FALSE)

# view result: PLS-DA plot
result_plsda$interacitve_plsda
static_plsda <- result_plsda$static_plsda
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/02.Dimensionality_reduction/PLS-DA/static_PLSDA_plot.png",
    width = 1200, height = 800, res = 150
)
print(static_plsda)
dev.off()

# view result: PLS-DA loading plot
result_plsda$interactive_loadingPlot
static_loadingPlot <- result_plsda$static_loadingPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/02.Dimensionality_reduction/PLS-DA/static_loading_plot.png",
    width = 2000, height = 2000, res = 150
)
print(static_loadingPlot)
dev.off()


#########################
# HIERARCHICAL CLUSTERING
# conduct hierarchical clustering by class
result_hcluster <- heatmap_clustering(
  de_se=deSp_se, char='class', distfun='pearson', 
  hclustfun='ward.D2', type='sig')

# result summary
summary(result_hcluster)

# view result: heatmap of significant lipid species
result_hcluster$interactive_heatmap
static_heatmap <- result_hcluster$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/03.Hierarchical_clustering/static_Lipid_species_hierarchical_clustering_heatmap_class_pearson_wardD2_sig.png",
    width = 2000, height = 2000, res = 150
)
print(static_heatmap)
dev.off()

# conduct hierarchical clustering
result_hcluster <- heatmap_clustering(
  de_se=deSp_se, char='class', distfun='spearman', 
  hclustfun='ward.D2', type='sig')

# view result: heatmap of significant lipid species
result_hcluster$interactive_heatmap
static_heatmap <- result_hcluster$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/03.Hierarchical_clustering/static_Lipid_species_hierarchical_clustering_heatmap_class_spearman_wardD2_sig.png",
    width = 2000, height = 2000, res = 150
)
print(static_heatmap)
dev.off()



# conduct hierarchical clustering by category
result_hcluster <- heatmap_clustering(
  de_se=deSp_se, char='Category', distfun='pearson', 
  hclustfun='ward.D2', type='sig')

# result summary
summary(result_hcluster)

# view result: heatmap of significant lipid species
result_hcluster$interactive_heatmap
static_heatmap <- result_hcluster$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/03.Hierarchical_clustering/static_Lipid_species_hierarchical_clustering_heatmap_category_pearson_wardD2_sig.png",
    width = 2000, height = 2000, res = 150
)
print(static_heatmap)
dev.off()


##########################
# CHARACTERISTICS ANALYSIS
# conduct characteristic analysis
result_char <- char_association(deSp_se, char='Category')

# result summary
summary(result_char)

# view result: bar chart
result_char$interactive_barPlot
static_barPlot <- result_char$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Barchart_of_significant_groups_function.png",
    width = 2000, height = 2000, res = 150
)
print(static_barPlot)
dev.off()

# view result: lollipop plot
result_char$interacitve_lollipop
static_lollipop <- result_char$static_lollipop
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Lollipop_chart_of_all_significant_groups_category.png",
    width = 2000, height = 2000, res = 150
)
print(static_lollipop)
dev.off()



# conduct characteristic analysis
result_char <- char_association(deSp_se, char='class')

# result summary
summary(result_char)

# view result: bar chart
result_char$interactive_barPlot
static_barPlot <- result_char$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Barchart_of_significant_classes.png",
    width = 800, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# view result: lollipop plot
result_char$interacitve_lollipop
static_lollipop <- result_char$static_lollipop
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Lollipop_chart_of_all_significant_classes.png",
    width = 800, height = 800, res = 150
)
print(static_lollipop)
dev.off()



###################################
# DOUBLE BOND-CHAIN LENGTH ANALYSIS
# conduct double bond-chain length analysis (without setting `char_feature`)
heatmap_all <- heatmap_chain_db(
  processed_se, char='class', char_feature=NULL, ref_group='healthy', 
  test='One-way ANOVA', significant='pval', p_cutoff=0.05, 
  FC_cutoff=1, transform='log10')

# result summary 
summary(heatmap_all)

# summary of total chain result
summary(heatmap_all$total_chain)

# view result: heatmap of total chain
heatmap_all$total_chain$static_heatmap
total_chain <- heatmap_all$total_chain$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Heatmap_of_total_chain_class.png",
    width = 2400, height = 1400, res = 150
)
print(total_chain)
dev.off()

write_tsv(heatmap_all$total_chain$table_heatmap, file = "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/total_chain_table_heatmap.tsv")
write_tsv(heatmap_all$total_chain$processed_abundance, file = "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/total_chain_processed_abundance.tsv")
write_tsv(heatmap_all$total_chain$transformed_abundance, file = "/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/total_chain_transformed_abundance.tsv")


# summary of each chain result
summary(heatmap_all$each_chain)

# view result: heatmap of each chain
heatmap_all$each_chain$static_heatmap
each_chain <- heatmap_all$each_chain$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Heatmap_of_each_chain_class.png",
    width = 2400, height = 800, res = 150
)
print(total_chain)
dev.off()

###########
# PC
# conduct double bond-chain length analysis (a specific `char_feature`)
heatmap_one <- heatmap_chain_db(
  processed_se, char='class', char_feature='PC', ref_group='healthy',
  test='One-way ANOVA', significant='padj', p_cutoff=0.05, 
  FC_cutoff=1, transform='log10')

# result summary 
summary(heatmap_one)

# summary of total chain result
summary(heatmap_one$total_chain)

# view result: heatmap of total chain
heatmap_one$total_chain$static_heatmap
total_chain <- heatmap_all$total_chain$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Heatmap_of_total_chain_class_PC.png",
    width = 2400, height = 1400, res = 150
)
print(total_chain)
dev.off()

# summary of each chain result
summary(heatmap_one$each_chain)

# view result: heatmap of each chain
heatmap_one$each_chain$static_heatmap
each_chain <- heatmap_all$each_chain$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Heatmap_of_each_chain_class_PC.png",
    width = 2400, height = 800, res = 150
)
print(total_chain)
dev.off()

# plot abundance box plot of "38:7"
# boxPlot_result <- boxPlot_feature_twoGroup(
#   heatmap_one$each_chain$chain_db_se, feature='38:7', 
#   ref_group='mild', test='t-test', transform='log10')
# 
# # result summary
# summary(boxPlot_result)
# 
# # view result: static box plot
# boxPlot_result$static_boxPlot

###########
# FA
# conduct double bond-chain length analysis (a specific `char_feature`)
heatmap_one <- heatmap_chain_db(
  processed_se, char='class', char_feature='FA', ref_group='mild',
  test='t-test', significant='pval', p_cutoff=0.05, 
  FC_cutoff=1, transform='log10')

# result summary 
summary(heatmap_one)

# summary of total chain result
summary(heatmap_one$total_chain)

# view result: heatmap of total chain
heatmap_one$total_chain$static_heatmap
total_chain <- heatmap_all$total_chain$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Heatmap_of_total_chain_class_FA.png",
    width = 2400, height = 1400, res = 150
)
print(total_chain)
dev.off()

# summary of each chain result
summary(heatmap_one$each_chain)

# view result: heatmap of each chain
heatmap_one$each_chain$static_heatmap
each_chain <- heatmap_all$each_chain$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Heatmap_of_each_chain_class_FA.png",
    width = 2400, height = 800, res = 150
)
print(total_chain)
dev.off()


###########
# LPC
# conduct double bond-chain length analysis (a specific `char_feature`)
heatmap_one <- heatmap_chain_db(
  processed_se, char='class', char_feature='LPC', ref_group='mild',
  test='t-test', significant='pval', p_cutoff=0.05, 
  FC_cutoff=1, transform='log10')

# result summary 
summary(heatmap_one)

# summary of total chain result
summary(heatmap_one$total_chain)

# view result: heatmap of total chain
heatmap_one$total_chain$static_heatmap
total_chain <- heatmap_all$total_chain$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Heatmap_of_total_chain_class_LPC.png",
    width = 2400, height = 1400, res = 150
)
print(total_chain)
dev.off()

# summary of each chain result
summary(heatmap_one$each_chain)

# view result: heatmap of each chain
heatmap_one$each_chain$static_heatmap
each_chain <- heatmap_all$each_chain$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Healthy_vs_sick_patients/D0/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Heatmap_of_each_chain_class_LPC.png",
    width = 2400, height = 800, res = 150
)
print(total_chain)
dev.off()
