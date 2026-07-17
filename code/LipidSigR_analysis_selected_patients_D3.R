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

group_info_twoGroup <- read_tsv("/Users/loictalignani/research/project/lipidomics/data/lipidsig_datasets/selected_patients/group_information_table_selected_patients_D3.tsv")
head(group_info_twoGroup, 5)

abundance_twoGroup <- read_tsv("/Users/loictalignani/research/project/lipidomics/data/lipidsig_datasets/selected_patients/MISSE_Global_Lipidomics_results_selected_patients_D3.tsv")
head(abundance_twoGroup)

# map lipid characteristics by rgoslin
library(rgoslin)
library(dplyr)

parse_lipid <- rgoslin::parseLipidNames(lipidNames=abundance_twoGroup$feature)

# filter lipid recognized by rgoslin
recognized_lipid <- parse_lipid$Original.Name[
  which(parse_lipid$Grammar != 'NOT_PARSEABLE')]
abundance <- abundance_twoGroup %>% 
  dplyr::filter(feature %in% recognized_lipid)
goslin_annotation <- parse_lipid %>% 
  dplyr::filter(Original.Name %in% recognized_lipid)

head(abundance[, 1:6], 5)

head(goslin_annotation[, 1:6], 5)

# Construct SE object
se <- as_summarized_experiment(
  abundance, goslin_annotation, group_info=group_info_twoGroup, 
  se_type='de_two', paired_sample=FALSE)

se@assays@data@listData

# data processing
processed_se <- data_process(
  se, exclude_missing=TRUE, exclude_missing_pct=50, 
  replace_na_method='min', replace_na_method_ref=0.5, 
  normalization='Percentage', transform='log10')

processed_se@metadata$processed_abund


#########################
# PROFILING
#########################

# conduct differential expression analysis of lipid species
deSp_se <- deSp_twoGroup(
  processed_se, ref_group='mild', test='t-test',
  significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')

# extract results in SE
res_list <- extract_summarized_experiment(deSp_se)
# summary of extract results
summary(res_list)

# DATA PROCESSING
# plotting
data_process_plots <- plot_data_process(se, processed_se)

# result summary
summary(data_process_plots)

# view box plot before/after data processing
data_process_plots$static_boxPlot_before


static_boxPlot_before <- data_process_plots$static_boxPlot_before
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/01.Preprocessing/00.Data_quality/BoxPlot_before_process.pdf")
print(data_process_plots$static_boxPlot_before)
dev.off()

data_process_plots$static_boxPlot_after

static_boxPlot_after <- data_process_plots$static_boxPlot_after
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/01.Preprocessing/00.Data_quality/Boxplot_after_process.pdf")
print(static_boxPlot_before)
dev.off()

data_process_plots$static_densityPlot_before

# view density plot before/after data processing
static_densityPlot_before <- data_process_plots$static_densityPlot_before
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/01.Preprocessing/00.Data_quality/Densityplot_before_process.pdf")
print(static_densityPlot_before)
dev.off()

data_process_plots$static_densityPlot_after

static_densityPlot_after <- data_process_plots$static_densityPlot_after  
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/01.Preprocessing/00.Data_quality/Densityplot_after_process.pdf")
print(static_densityPlot_after)
dev.off()


# CROSS-SAMPLE VARIABILITY
## conduct profiling
result <- cross_sample_variability(se)

## result summary
summary(result)

## view result: histogram of lipid numbers
result$static_lipid_number_barPlot
static_lipid_number_barPlot <- result$static_lipid_number_barPlot
png(
  "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/01.Preprocessing/01.Cross-sample_variability/Expressed_lipid_numbers.png",
  width = 1600, height = 1000, res = 150
)
print(static_lipid_number_barPlot)
dev.off()

# view result: histogram of the total amount of lipid in each sample.
result$static_lipid_amount_barPlot
static_lipid_amount_barPlot <- result$static_lipid_amount_barPlot
png(
  "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/01.Preprocessing/01.Cross-sample_variability/Histogram_lipid_amount.png",
  width = 1400, height = 1000, res = 150
)
print(static_lipid_amount_barPlot)
dev.off()

# view result: density plot of the underlying probability distribution
result$static_lipid_distribution
static_lipid_distribution <- result$static_lipid_distribution
png(
  "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/01.Preprocessing/01.Cross-sample_variability/lipid_abundance_distribution.png",
  width = 1400, height = 1000, res = 150
)
print(static_lipid_distribution)
dev.off()


# DIMENSIONALITY REDUCTION

## PCA
# conduct PCA
result_pca <- dr_pca(
  processed_se, scaling=TRUE, centering=TRUE, clustering='kmeans', 
  cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, 
  eps=NULL, minPts=NULL, feature_contrib_pc=c(1,2), plot_topN=10)

# result summary
summary(result_pca)

# view result: PCA plot
result_pca$interactive_pca
static_pca <- result_pca$static_pca
png(
  "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/02.Results_D3/01.Dimensionality_reduction/PCA/static_PCA_plot.png",
  width = 1400, height = 1000, res = 150
)
print(static_pca)
dev.off()

# view result: scree plot of top 10 principle components
result_pca$static_screePlot
static_screePlot <- result_pca$static_screePlot
png(
  "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/02.Results_D3/01.Dimensionality_reduction/PCA/Static_explained_variance.png",
  width = 1400, height = 1000, res = 150
)
print(static_screePlot)
dev.off()

# view result: correlation circle plot of PCA variables
result_pca$static_feature_contribution
static_feature_contribution <- result_pca$static_feature_contribution
png(
  "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/02.Results_D3/01.Dimensionality_reduction/PCA/static_feature_contribution.png",
  width = 1400, height = 1000, res = 150
)
print(static_feature_contribution)
dev.off()

# view result: Correlation of contribution of top 10 features
result_pca$static_variablePlot
static_variablePlot <- result_pca$static_variablePlot
png(
  "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/02.Results_D3/01.Dimensionality_reduction/PCA/corr_circle_plot.png",
  width = 1400, height = 1000, res = 150
)
print(static_variablePlot)
dev.off()

# save results
write_tsv(result_pca$table_pca_contribution, file = "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/02.Results_D3/01.Dimensionality_reduction/PCA/table_pca_contribution.tsv")
write_tsv(result_pca$pca_rotated_data, file = "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/02.Results_D3/01.Dimensionality_reduction/PCA/table_pca_rotated_data.tsv")
# T-SNE


# conduct t-SNE
result_tsne <- dr_tsne(
  processed_se, pca=TRUE, perplexity=5, max_iter=500, clustering='kmeans',
  cluster_num=2, kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, 
  eps=NULL, minPts=NULL)

# result summary
summary(result_tsne)

# view result: t-SNE plot
result_tsne$interactive_tsne
static_tsne <- result_tsne$static_tsne
png(
  "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/02.Results_D3/01.Dimensionality_reduction/t-SNE/static_t-SNE_plot.png",
  width = 1400, height = 1000, res = 150
)
print(static_tsne)
dev.off()

# save results
write_tsv(result_tsne$tsne_result, file = "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/02.Results_D3/01.Dimensionality_reduction/t-SNE/table_tsne_results.tsv")


# UMAP
# conduct UMAP
result_umap <- dr_umap(
  processed_se, n_neighbors=15, scaling=TRUE, umap_metric='euclidean',
  clustering='kmeans', cluster_num=2, kmedoids_metric=NULL,
  distfun=NULL, hclustfun=NULL, eps=NULL, minPts=NULL)

# result summary
summary(result_umap)

# view result: UMAP plot
result_umap$static_umap
static_umap <- result_umap$static_umap
png(
  "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/02.Results_D3/01.Dimensionality_reduction/UMAP/static_UMAP_plot.png",
  width = 1400, height = 1000, res = 150
)
print(static_umap)
dev.off()

# save umap results
write_tsv(result_umap$umap_result, file = "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/02.Results_D3/01.Dimensionality_reduction/UMAP/UMAP_results.tsv")

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
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/02.Results_D3/02.Correlation_Heatmap/by_samples/Correlation_Heatmap_by_samples_spearman_Eucl_WD2.png",
  width = 1400, height = 1000, res = 150
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
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/01.Profiling/02.Results_D3/02.Correlation_Heatmap/by_category/Correlation_Heatmap_by_lipid_category_spearman_Eucl_WD2.png",
    width = 1500, height = 1500, res = 150
)
print(static_heatmap)
dev.off()

# correlation calculation by class and by sample
result_heatmap <- heatmap_correlation(
  processed_se, char="class", transform='log10', correlation='spearman', 
  distfun='euclidean', hclustfun='ward.D2', type='sample')

result_heatmap$static_heatmap


# LIPID CHARACTERISTICS PER SAMPLE
# calculate lipid expression of selected characteristic
result_lipid <- lipid_profiling(processed_se, char="class")

# result summary
summary(result_lipid)

# view result: bar plot of class per sample
result_lipid$interactive_char_barPlot
static_char_barPlot <- result_lipid$static_char_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_expression/Lipid_class_per_sample.png",
    width = 1000, height = 1000, res = 150
)
print(static_char_barPlot)
dev.off()

# view result: stacked horizontal bar chart
result_lipid$interactive_lipid_composition
static_lipid_composition <- result_lipid$static_lipid_composition
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_expression/Lipid_composition_per_sample.png",
    width = 1000, height = 1000, res = 150
)
print(static_lipid_composition)
dev.off()


 # calculate lipid expression of selected characteristic: Total.C
totalC_profile <- lipid_profiling(processed_se, char = "Total.C")
totalC_profile$interactive_char_barPlot
static_char_barPlot <- totalC_profile$static_char_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_expression/Total.C_per_sample.png",
    width = 1000, height = 1000, res = 150
)
print(static_char_barPlot)
dev.off()


# calculate lipid expression of selected characteristic: Total.DB
totalDB_profile <- lipid_profiling(processed_se, char = "Total.DB")
totalDB_profile$interactive_char_barPlot
static_char_barPlot <- totalDB_profile$static_char_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_expression/Total.DB_per_sample.png",
    width = 1000, height = 1000, res = 150
)
print(static_char_barPlot)
dev.off()

# calculate lipid expression of selected characteristic: Total.OH
totalOH_profile <- lipid_profiling(processed_se, char = "Total.OH")
totalOH_profile$interactive_char_barPlot
static_char_barPlot <- totalOH_profile$static_char_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_expression/Total.OH_per_sample.png",
    width = 1000, height = 1000, res = 150
)
print(static_char_barPlot)
dev.off()



# calculate lipid expression of selected characteristic: Total.FA
totalFA_profile <- lipid_profiling(processed_se, char = "Total.FA")
totalFA_profile$interactive_char_barPlot
static_char_barPlot <- totalFA_profile$static_char_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_expression/Total.FA_per_sample.png",
    width = 3000, height = 3000, res = 150
)
print(static_char_barPlot)
dev.off()


##################################################
# DIFFERENTIAL EXPRESSION ANALYSIS
##################################################

# STEP 1: CONDUCT DIFFERENTIAL EXPRESSION
# conduct differential expression analysis of lipid species
deSp_se <- deSp_twoGroup(
  processed_se, ref_group='mild', test='t-test',
  significant='pval', p_cutoff=0.05, FC_cutoff=1, transform='log10')


# plot differential expression analysis result
deSp_plot <- plot_deSp_twoGroup(deSp_se)

# result summary
summary(deSp_plot)

# view result: lollipop chart
deSp_plot$static_de_lipid
static_de_lipid <- deSp_plot$static_de_lipid
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lollipop_chart/lollipop_chart.png",
    width = 1000, height = 1000, res = 150
)
print(static_de_lipid)
dev.off()

# view result: MA plot
deSp_plot$interactive_maPlot
static_maPlot <- deSp_plot$static_maPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/MA_plot/MA_plot.png",
    width = 1000, height = 1000, res = 150
)
print(static_maPlot)
dev.off()

# view result: Volcano plot
deSp_plot$interactive_volcanoPlot
static_volcanoPlot <- deSp_plot$static_volcanoPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Volcano_plot/Volcano_plot.png",
    width = 1000, height = 1000, res = 150
)
print(static_volcanoPlot)
dev.off()

write_tsv(deSp_plot$table_de_lipid, file = "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lollipop_chart/table_diffexp.tsv")
write_tsv(deSp_plot$table_ma_volcano, file = "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/MA_plot/table_MA_volcano.tsv")



# LIPID ABUNDANCE
# plot abundance box plot of 'CE 18:2' ---> this function doesn't work. I used the web site
# boxPlot_result <- boxPlot_feature_twoGroup(
#   processed_se, feature='CE 18:2', ref_group='mild', test='t-test',
#   transform='log10')
# 
# # result summary
# summary(boxPlot_result)
# 
# # view result: static box plot
# boxPlot_result$static_boxPlot



# LIPID CHARACTERISTICS DIFFERENTIAL EXPRESSION ANALYSIS
# two way anova
twoWayAnova_table <- char_2wayAnova(
  processed_se, ratio_transform='log2', char_transform='log10')

# view result table to select an available characteristic as char for deChar_twoGroup() function
head(twoWayAnova_table[, 1:4], 5) 

list_lipid_char(processed_se)

# conduct differential expression of lipid characteristics for Fatty acid properties "Total.C"
deChar_se <- deChar_twoGroup(
  processed_se, char="Total.C", ref_group="mild", test='t-test', 
  significant="pval", p_cutoff=0.05, FC_cutoff=1, transform='log10')

# plot differential expression analysis results
deChar_plot <- plot_deChar_twoGroup(deChar_se)

# result summary
summary(deChar_plot)

# view result: bar plot of selected `Total.C`
deChar_plot$interactive_barPlot
static_barPlot <- deChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/TotalC/Barplot_of_characteristics_analysis_TotalC.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# conduct differential expression of lipid characteristics for Fatty acid properties "Total.DB"
deChar_se <- deChar_twoGroup(
  processed_se, char="Total.DB", ref_group="mild", test='t-test', 
  significant="pval", p_cutoff=0.05, FC_cutoff=1, transform='log10')

# plot differential expression analysis results
deChar_plot <- plot_deChar_twoGroup(deChar_se)

# result summary
summary(deChar_plot)

# view result: bar plot of selected `Total.DB`
deChar_plot$interactive_barPlot
static_barPlot <- deChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/TotalDB/Barplot_of_characteristics_analysis_TotalDB.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# conduct differential expression of lipid characteristics for Fatty acid properties "Total.FA"
deChar_se <- deChar_twoGroup(
  processed_se, char="Total.FA", ref_group="mild", test='t-test', 
  significant="pval", p_cutoff=0.05, FC_cutoff=1, transform='log10')

# plot differential expression analysis results
deChar_plot <- plot_deChar_twoGroup(deChar_se)

# result summary
summary(deChar_plot)

# view result: bar plot of selected `Total.FA`
deChar_plot$interactive_barPlot
static_barPlot <- deChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/TotalFA/Barplot_of_characteristics_analysis_TotalFA.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# conduct differential expression of lipid characteristics for Fatty acid properties "Total.OH"
deChar_se <- deChar_twoGroup(
  processed_se, char="Total.OH", ref_group="mild", test='t-test', 
  significant="pval", p_cutoff=0.05, FC_cutoff=1, transform='log10')

# plot differential expression analysis results
deChar_plot <- plot_deChar_twoGroup(deChar_se)

# result summary
summary(deChar_plot)

# view result: bar plot of selected `Total.OH`
deChar_plot$interactive_barPlot
static_barPlot <- deChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/TotalOH/Barplot_of_characteristics_analysis_TotalOH.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# conduct differential expression of lipid characteristics for "class"
deChar_se <- deChar_twoGroup(
  processed_se, char="class", ref_group="mild", test='t-test', 
  significant="pval", p_cutoff=0.05, FC_cutoff=1, transform='log10')

# plot differential expression analysis results
deChar_plot <- plot_deChar_twoGroup(deChar_se)

# result summary
summary(deChar_plot)

# view result: bar plot of selected `char`
deChar_plot$interactive_barPlot
static_barPlot <- deChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Class/Barplot_of_characteristics_analysis_class.png",
    width = 1200, height = 800, res = 150
    )
print(static_barPlot)
dev.off()


# conduct differential expression of lipid characteristics for "Category"
deChar_se <- deChar_twoGroup(
  processed_se, char="Category", ref_group="mild", test='t-test', 
  significant="pval", p_cutoff=0.05, FC_cutoff=1, transform='log10')

# plot differential expression analysis results
deChar_plot <- plot_deChar_twoGroup(deChar_se)

# result summary
summary(deChar_plot)

# view result: bar plot of selected `char`
deChar_plot$interactive_barPlot
static_barPlot <- deChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Category/Barplot_of_characteristics_analysis_category.png",
    width = 1200, height = 800, res = 150
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
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_LPCO-.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="TG")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_TG.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PC")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_PC.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="CE")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_CE.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="DG")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_DG.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="LPC")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_LPC.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="FA")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_FA.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="CAR")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_CAR.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PI")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_PI.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="Cer")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_Cer.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PC O-")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalC/Total_length_per_PC O-.png",
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
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_LPCO-.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="TG")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_TG.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PC")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_PC.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="CE")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_CE.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="DG")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_DG.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="LPC")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_LPC.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="FA")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_FA.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="CAR")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_CAR.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PI")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_PI.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="Cer")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_Cer.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PC O-")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalDB/Total_DB_per_PC O-.png",
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
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_LPCO-.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="TG")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_TG.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PC")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_PC.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="CE")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_CE.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="DG")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_DG.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="LPC")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_LPC.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="FA")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_FA.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="CAR")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_CAR.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PI")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_PI.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="Cer")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_Cer.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="PC O-")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Species_based_aggregation/TotalOH/Total_OH_per_PC O-.png",
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
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Fatty_acids_based_aggregation/per_category/Glycerophospholipids.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="Sterol lipids [ST]")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Fatty_acids_based_aggregation/per_category/Sterol_lipids.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# visualize subgroup differential expression of lipid characteristics 
subChar_plot <- plot_subChar_twoGroup(subChar_se, subChar_feature="Fatty acyls [FA]")

# view result: bar plot of `subChar_feature`
subChar_plot$static_barPlot
static_barPlot <- subChar_plot$static_barPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step1/01.DiffExp/Lipid_characteristics_DiffExp/Fatty_acids_based_aggregation/per_category/Fatty_acyls.png",
    width = 1200, height = 800, res = 150
)
print(static_barPlot)
dev.off()



#####################################
# STEP 2: DIMENSION REDUCTION: PLS-DA
# conduct PLSDA
result_plsda <- dr_plsda(
  deSp_se, ncomp=2, scaling=TRUE, clustering='group_info', cluster_num=2, 
  kmedoids_metric=NULL, distfun=NULL, hclustfun=NULL, eps=NULL, minPts=NULL)

# result summary
summary(result_plsda)
write.csv(result_plsda$plsda_result, file = "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/02.Dimensionality_reduction/PLS-DA/Table_of_sample_variate_D3.csv", row.names = FALSE)
write.csv(result_plsda$table_plsda_loading, file = "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/02.Dimensionality_reduction/PLS-DA/Table_of_sample_loading_D3.csv", row.names = FALSE)

# view result: PLS-DA plot
result_plsda$interacitve_plsda
static_plsda <- result_plsda$static_plsda
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/02.Dimensionality_reduction/PLS-DA/static_PLSDA_plot.png",
    width = 1200, height = 800, res = 150
)
print(static_plsda)
dev.off()

# view result: PLS-DA loading plot
result_plsda$interactive_loadingPlot
static_loadingPlot <- result_plsda$static_loadingPlot
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/02.Dimensionality_reduction/PLS-DA/static_loading_plot.png",
    width = 800, height = 800, res = 150
)
print(static_loadingPlot)
dev.off()


#########################
# HIERARCHICAL CLUSTERING
# conduct hierarchical clustering
result_hcluster <- heatmap_clustering(
  de_se=deSp_se, char='class', distfun='pearson', 
  hclustfun='ward.D2', type='sig')

# result summary
summary(result_hcluster)

# view result: heatmap of significant lipid species
result_hcluster$interactive_heatmap
static_heatmap <- result_hcluster$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/03.Hierarchical_clustering/static_Lipid_species_hierarchical_clustering_heatmap_class_pearson_wardD2_sig.png",
    width = 800, height = 800, res = 150
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
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Barchart_of_significant_groups_function.png",
    width = 800, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# view result: lollipop plot
result_char$interacitve_lollipop
static_lollipop <- result_char$static_lollipop
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Lollipop_chart_of_all_significant_groups_category.png",
    width = 800, height = 800, res = 150
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
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Barchart_of_significant_classes.png",
    width = 800, height = 800, res = 150
)
print(static_barPlot)
dev.off()

# view result: lollipop plot
result_char$interacitve_lollipop
static_lollipop <- result_char$static_lollipop
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Lollipop_chart_of_all_significant_classes.png",
    width = 800, height = 800, res = 150
)
print(static_lollipop)
dev.off()



###################################
# DOUBLE BOND-CHAIN LENGTH ANALYSIS
# conduct double bond-chain length analysis (without setting `char_feature`)
heatmap_all <- heatmap_chain_db(
  processed_se, char='class', char_feature=NULL, ref_group='mild', 
  test='t-test', significant='pval', p_cutoff=0.05, 
  FC_cutoff=1, transform='log10')

# result summary 
summary(heatmap_all)

# summary of total chain result
summary(heatmap_all$total_chain)

# view result: heatmap of total chain
heatmap_all$total_chain$static_heatmap
total_chain <- heatmap_all$total_chain$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Heatmap_of_total_chain_class.png",
    width = 2400, height = 1400, res = 150
)
print(total_chain)
dev.off()

write_tsv(heatmap_all$total_chain$table_heatmap, file = "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/total_chain_table_heatmap.tsv")
write_tsv(heatmap_all$total_chain$processed_abundance, file = "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/total_chain_processed_abundance.tsv")
write_tsv(heatmap_all$total_chain$transformed_abundance, file = "/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/total_chain_transformed_abundance.tsv")


# summary of each chain result
summary(heatmap_all$each_chain)

# view result: heatmap of each chain
heatmap_all$each_chain$static_heatmap
each_chain <- heatmap_all$each_chain$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Heatmap_of_each_chain_class.png",
    width = 2400, height = 800, res = 150
)
print(total_chain)
dev.off()

###########
# PC
# conduct double bond-chain length analysis (a specific `char_feature`)
heatmap_one <- heatmap_chain_db(
  processed_se, char='class', char_feature='PC', ref_group='mild',
  test='t-test', significant='pval', p_cutoff=0.05, 
  FC_cutoff=1, transform='log10')

# result summary 
summary(heatmap_one)

# summary of total chain result
summary(heatmap_one$total_chain)

# view result: heatmap of total chain
heatmap_one$total_chain$static_heatmap
total_chain <- heatmap_all$total_chain$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Heatmap_of_total_chain_class_PC.png",
    width = 2400, height = 1400, res = 150
)
print(total_chain)
dev.off()

# summary of each chain result
summary(heatmap_one$each_chain)

# view result: heatmap of each chain
heatmap_one$each_chain$static_heatmap
each_chain <- heatmap_all$each_chain$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Heatmap_of_each_chain_class_PC.png",
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
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Heatmap_of_total_chain_class_FA.png",
    width = 2400, height = 1400, res = 150
)
print(total_chain)
dev.off()

# summary of each chain result
summary(heatmap_one$each_chain)

# view result: heatmap of each chain
heatmap_one$each_chain$static_heatmap
each_chain <- heatmap_all$each_chain$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Heatmap_of_each_chain_class_FA.png",
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
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Heatmap_of_total_chain_class_LPC.png",
    width = 2400, height = 1400, res = 150
)
print(total_chain)
dev.off()

# summary of each chain result
summary(heatmap_one$each_chain)

# view result: heatmap of each chain
heatmap_one$each_chain$static_heatmap
each_chain <- heatmap_all$each_chain$static_heatmap
png("/Users/loictalignani/research/project/lipidomics/analysis/Selected_patients/D3/02.DiffExp/02.DiffExp_analysis/Lipid_species_analysis/Step2/04.Characteristics_association/Heatmap_of_each_chain_class_LPC.png",
    width = 2400, height = 800, res = 150
)
print(total_chain)
dev.off()
