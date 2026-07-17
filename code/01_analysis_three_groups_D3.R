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

# Set working directory
setwd("E:/PBMC_lipidomics")

# Define base directory for outputs
base_dir <- "/Users/loictalignani/research/project/lipidomics/analysis/Three_groups/D3"

cat("\n========================================\n")
cat("THREE-GROUP LIPIDOMICS ANALYSIS\n")
cat("Groups: Healthy, Mild, Severe\n")
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
dir.create(file.path(base_dir, "02.DiffExp/01.ANOVA"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "02.DiffExp/02.Visualizations"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "02.DiffExp/03.Individual_lipid_boxplots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "03.Enrichment/LSEA"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "03.Enrichment/ORA"), recursive = TRUE, showWarnings = FALSE)

# Load group information
cat("Loading group information...\n")
group_info <- read_tsv("data/lipidsig_datasets/healthy_vs_sick_patients/group_information_table_healthy_vs_sick_patients_D3.tsv")
cat("  Groups found:", unique(group_info$group), "\n")
cat("  Sample count:", nrow(group_info), "\n")

# Load abundance data
cat("Loading abundance data...\n")
abundance <- read_tsv("data/lipidsig_datasets/healthy_vs_sick_patients/healthy_sick_lipidomics_D3.tsv")
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
  se_type = "de_multiple", # Multi-group comparison
  paired_sample = NULL
)

cat("  ✓ SummarizedExperiment created\n")

# Data processing
cat("Processing data (filtering, normalization, transformation)...\n")
processed_se <- data_process(
  se,
  exclude_missing = TRUE, # Remove features with too many missing values
  exclude_missing_pct = 70, # Keep features present in ≥70% samples
  replace_na_method = "min", # Replace NA with minimum value
  replace_na_method_ref = 0.5, # Factor for minimum value
  normalization = "Percentage", # Normalization method
  transform = "log10" # Log10 transformation
)

cat("  ✓ Data processing complete\n")


##########################################################################################
# 2. PROFILING - DATA QUALITY
##########################################################################################

cat("\n=== 2. PROFILING - DATA QUALITY ===\n")

# Plot data processing effects
cat("Generating data quality plots...\n")
data_process_plots <- plot_data_process(se, processed_se)

png(file.path(base_dir, "01.Profiling/00.Data_quality/BoxPlot_before_process.png"),
  width = 1400, height = 1000, res = 150
)
print(data_process_plots$static_boxPlot_before)
dev.off()

png(file.path(base_dir, "01.Profiling/00.Data_quality/BoxPlot_after_process.png"),
  width = 1400, height = 1000, res = 150
)
print(data_process_plots$static_boxPlot_after)
dev.off()

png(file.path(base_dir, "01.Profiling/00.Data_quality/DensityPlot_before_process.png"),
  width = 1400, height = 1000, res = 150
)
print(data_process_plots$static_densityPlot_before)
dev.off()

png(file.path(base_dir, "01.Profiling/00.Data_quality/DensityPlot_after_process.png"),
  width = 1400, height = 1000, res = 150
)
print(data_process_plots$static_densityPlot_after)
dev.off()

cat("  ✓ Data quality plots saved\n")


# Cross-sample variability
cat("Analyzing cross-sample variability...\n")
cross_var_result <- cross_sample_variability(se)

png(file.path(base_dir, "01.Profiling/01.Cross-sample_variability/Expressed_lipid_numbers.png"),
  width = 1600, height = 1000, res = 150
)
print(cross_var_result$static_lipid_number_barPlot)
dev.off()

png(file.path(base_dir, "01.Profiling/01.Cross-sample_variability/Histogram_lipid_amount.png"),
  width = 1400, height = 1000, res = 150
)
print(cross_var_result$static_lipid_amount_barPlot)
dev.off()

png(file.path(base_dir, "01.Profiling/01.Cross-sample_variability/Lipid_abundance_distribution.png"),
  width = 1400, height = 1000, res = 150
)
print(cross_var_result$static_lipid_distribution)
dev.off()

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

png(file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/PCA_plot.png"),
  width = 1400, height = 1000, res = 150
)
print(result_pca$static_pca)
dev.off()

# Save interactive PCA plot as RDS
if (!is.null(result_pca$interactive_pca)) {
  saveRDS(
    result_pca$interactive_pca,
    file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/PCA_plot_interactive.rds")
  )
}

png(file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Scree_plot.png"),
  width = 1400, height = 1000, res = 150
)
print(result_pca$static_screePlot)
dev.off()

# Save interactive scree plot as RDS
if (!is.null(result_pca$interactive_screePlot)) {
  saveRDS(
    result_pca$interactive_screePlot,
    file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Scree_plot_interactive.rds")
  )
}

png(file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Feature_contribution.png"),
  width = 1400, height = 1000, res = 150
)
print(result_pca$static_feature_contribution)
dev.off()

# Save interactive feature contribution as RDS
if (!is.null(result_pca$interactive_feature_contribution)) {
  saveRDS(
    result_pca$interactive_feature_contribution,
    file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Feature_contribution_interactive.rds")
  )
}

png(file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Variable_correlation.png"),
  width = 1400, height = 1000, res = 150
)
print(result_pca$static_variablePlot)
dev.off()

# Save interactive variable correlation as RDS
if (!is.null(result_pca$interactive_variablePlot)) {
  saveRDS(
    result_pca$interactive_variablePlot,
    file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/PCA/Variable_correlation_interactive.rds")
  )
}

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

png(file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/t-SNE/t-SNE_plot.png"),
  width = 1400, height = 1000, res = 150
)
print(result_tsne$static_tsne)
dev.off()

# Save interactive t-SNE plot as RDS
if (!is.null(result_tsne$interactive_tsne)) {
  saveRDS(
    result_tsne$interactive_tsne,
    file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/t-SNE/t-SNE_plot_interactive.rds")
  )
}

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

png(file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/UMAP/UMAP_plot.png"),
  width = 1400, height = 1000, res = 150
)
print(result_umap$static_umap)
dev.off()

# Save interactive UMAP plot as RDS
if (!is.null(result_umap$interactive_umap)) {
  saveRDS(
    result_umap$interactive_umap,
    file.path(base_dir, "01.Profiling/02.Dimensionality_reduction/UMAP/UMAP_plot_interactive.rds")
  )
}

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

png(file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_samples/Heatmap_by_samples.png"),
  width = 2000, height = 2000, res = 150
)
print(result_heatmap_sample$static_heatmap)
dev.off()

# Save interactive heatmap as RDS
if (!is.null(result_heatmap_sample$interactive_heatmap)) {
  saveRDS(
    result_heatmap_sample$interactive_heatmap,
    file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_samples/Heatmap_by_samples_interactive.rds")
  )
}

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

png(file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_category/Heatmap_by_category.png"),
  width = 1500, height = 1500, res = 150
)
print(result_heatmap_category$static_heatmap)
dev.off()

# Save interactive heatmap as RDS
if (!is.null(result_heatmap_category$interactive_heatmap)) {
  saveRDS(
    result_heatmap_category$interactive_heatmap,
    file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_category/Heatmap_by_category_interactive.rds")
  )
}

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

png(file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_class/Heatmap_by_class.png"),
  width = 1500, height = 1500, res = 150
)
print(result_heatmap_class$static_heatmap)
dev.off()

# Save interactive heatmap as RDS
if (!is.null(result_heatmap_class$interactive_heatmap)) {
  saveRDS(
    result_heatmap_class$interactive_heatmap,
    file.path(base_dir, "01.Profiling/03.Correlation_Heatmap/by_class/Heatmap_by_class_interactive.rds")
  )
}

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

png(file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Lipid_class_per_sample.png"),
  width = 1000, height = 1000, res = 150
)
print(result_lipid_class$static_char_barPlot)
dev.off()

png(file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Lipid_composition_per_sample.png"),
  width = 1000, height = 1000, res = 150
)
print(result_lipid_class$static_lipid_composition)
dev.off()


# Lipid profiling by Total Carbon
cat("Analyzing lipid composition by Total.C...\n")
result_totalC <- lipid_profiling(processed_se, char = "Total.C")

png(file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Total_C_per_sample.png"),
  width = 1000, height = 1000, res = 150
)
print(result_totalC$static_char_barPlot)
dev.off()


# Lipid profiling by Total Double Bonds
cat("Analyzing lipid composition by Total.DB...\n")
result_totalDB <- lipid_profiling(processed_se, char = "Total.DB")

png(file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Total_DB_per_sample.png"),
  width = 1000, height = 1000, res = 150
)
print(result_totalDB$static_char_barPlot)
dev.off()


# Lipid profiling by Category
cat("Analyzing lipid composition by Category...\n")
result_category <- lipid_profiling(processed_se, char = "Category")

png(file.path(base_dir, "01.Profiling/04.Lipid_characteristics/Category_per_sample.png"),
  width = 1000, height = 1000, res = 150
)
print(result_category$static_char_barPlot)
dev.off()

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

  # Count significant lipids
  n_sig <- sum(desp_df$adj.P.Val < 0.05, na.rm = TRUE)
  cat("  Total lipids analyzed:", nrow(desp_df), "\n")
  cat("  Significant lipids (adj.P.Val < 0.05):", n_sig, "\n")

  # Save all results
  write_tsv(desp_df,
    file = file.path(base_dir, "02.DiffExp/01.ANOVA/ANOVA_all_results.tsv")
  )

  # Save significant results separately
  if (n_sig > 0) {
    sig_df <- desp_df %>% filter(adj.P.Val < 0.05)
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
png(file.path(base_dir, "02.DiffExp/02.Visualizations/Lollipop_chart.png"),
  width = 1000, height = 1000, res = 150
)
print(deSp_plot$static_de_lipid)
dev.off()

cat("  ✓ ANOVA analysis complete\n")


# Pairwise comparisons within ANOVA (if available in LipidSigR)
# This provides detailed comparisons between each pair of groups
cat("Generating pairwise comparison plots...\n")

# Generate heatmap of significantly different lipids
if (!is.null(deSp_plot$static_heatmap)) {
  png(file.path(base_dir, "02.DiffExp/02.Visualizations/Heatmap_significant_lipids.png"),
    width = 1500, height = 1200, res = 150
  )
  print(deSp_plot$static_heatmap)
  dev.off()

  # Save interactive heatmap as RDS
  if (!is.null(deSp_plot$interactive_heatmap)) {
    saveRDS(
      deSp_plot$interactive_heatmap,
      file.path(base_dir, "02.DiffExp/02.Visualizations/Heatmap_significant_lipids_interactive.rds")
    )
  }
}

cat("  ✓ Visualization complete\n")


# Individual lipid species boxplots
cat("\n--- Creating individual lipid abundance boxplots ---\n")

# Get significant lipids from ANOVA results
if (exists("desp_df") && !is.null(desp_df) && nrow(desp_df) > 0) {
  # Check which p-value column is available
  # Prioritize adj.P.Val to match the earlier counting
  pval_col <- if ("adj.P.Val" %in% colnames(desp_df)) {
    "adj.P.Val"
  } else if ("padj" %in% colnames(desp_df)) {
    "padj"
  } else if ("adj.pval" %in% colnames(desp_df)) {
    "adj.pval"
  } else if ("P.Value" %in% colnames(desp_df)) {
    "P.Value"
  } else {
    NULL
  }

  if (!is.null(pval_col)) {
    cat("  Using p-value column:", pval_col, "\n")

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

      # Extract abundance data and group information
      abundance_data <- as.data.frame(SummarizedExperiment::assay(processed_se))
      group_data <- as.data.frame(SummarizedExperiment::colData(processed_se))

      for (i in 1:max_plots) {
        # Get lipid name - check if column exists
        # Ensure we extract a single character string
        lipid_name <- if ("feature" %in% colnames(sig_lipids)) {
          as.character(sig_lipids$feature[i])
        } else if ("lipid" %in% colnames(sig_lipids)) {
          as.character(sig_lipids$lipid[i])
        } else {
          as.character(rownames(sig_lipids)[i])
        }

        # Ensure lipid_name is a single value
        if (length(lipid_name) != 1 || is.na(lipid_name)) {
          cat("    Warning: Invalid lipid name at position", i, "\n")
          next
        }

        # Create safe filename (replace special characters)
        safe_name <- gsub("[^A-Za-z0-9_-]", "_", lipid_name)

        tryCatch(
          {
            # Get p-value for this lipid
            pval <- sig_lipids[[pval_col]][i]

            # Check if lipid exists in abundance data
            if (!lipid_name %in% rownames(abundance_data)) {
              cat("    Warning: Lipid", lipid_name, "not found in abundance data\n")
              next
            }

            # Extract abundance values for this lipid
            lipid_values <- as.numeric(abundance_data[lipid_name, ])

            # Create data frame for plotting
            plot_data <- data.frame(
              Sample = colnames(abundance_data),
              Abundance = lipid_values,
              Group = group_data$group
            )

            # Create boxplot with ggplot2
            library(ggplot2)
            library(ggsignif)

            # Perform pairwise t-tests for significance bars
            pairwise_comparisons <- list(
              c("healthy", "mild"),
              c("healthy", "severe"),
              c("mild", "severe")
            )

            # Calculate p-values for each comparison
            pairwise_pvals <- sapply(pairwise_comparisons, function(comp) {
              group1_vals <- plot_data$Abundance[plot_data$Group == comp[1]]
              group2_vals <- plot_data$Abundance[plot_data$Group == comp[2]]
              if (length(group1_vals) > 0 && length(group2_vals) > 0) {
                t_test <- t.test(group1_vals, group2_vals)
                return(t_test$p.value)
              } else {
                return(1)
              }
            })

            # Create significance labels (stars)
            sig_labels <- sapply(pairwise_pvals, function(p) {
              if (p < 0.0001) return("****")
              else if (p < 0.001) return("***")
              else if (p < 0.01) return("**")
              else if (p < 0.05) return("*")
              else return("ns")
            })

            p <- ggplot(plot_data, aes(x = Group, y = Abundance, fill = Group)) +
              geom_boxplot(outlier.shape = NA, alpha = 0.7) +
              geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
              geom_signif(
                comparisons = pairwise_comparisons,
                annotations = sig_labels,
                map_signif_level = FALSE,
                textsize = 4,
                vjust = 0.5,
                step_increase = 0.1
              ) +
              labs(
                title = lipid_name,
                subtitle = sprintf("ANOVA adj.p-val = %.2e", pval),
                x = "Group",
                y = "Log10 Abundance"
              ) +
              theme_bw() +
              theme(
                plot.title = element_text(face = "bold", size = 14),
                plot.subtitle = element_text(size = 10),
                axis.title = element_text(size = 12),
                axis.text = element_text(size = 10),
                legend.position = "none"
              ) +
              scale_fill_manual(values = c(
                "healthy" = "#4DBBD5FF",
                "mild" = "#E64B35FF",
                "severe" = "#00A087FF"
              ))

            # Save the plot
            png(
              file.path(
                base_dir,
                "02.DiffExp/03.Individual_lipid_boxplots",
                paste0(sprintf("%03d", i), "_", safe_name, ".png")
              ),
              width = 800, height = 600, res = 120
            )
            print(p)
            dev.off()

            # Progress indicator every 10 lipids
            if (i %% 10 == 0) {
              cat("    Progress:", i, "/", max_plots, "\n")
            }
          },
          error = function(e) {
            cat(
              "    Warning: Could not create plot for", lipid_name, "-",
              conditionMessage(e), "\n"
            )
          }
        )
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
# 7. LIPID SET ENRICHMENT ANALYSIS (LSEA)
##########################################################################################

cat("\n=== 7. LIPID SET ENRICHMENT ANALYSIS (LSEA) ===\n")

# LSEA uses the differential expression results stored in deSp_se
# It ranks lipids and tests for enrichment in different lipid sets

# Perform LSEA using LipidSigR
cat("Running LSEA by all characteristics...\n")

# LSEA with all characteristics (comprehensive)
lsea_all_result <- enrichment_lsea(
  deSp_se,
  char = NULL, # NULL = all characteristics
  rank_by = "statistic", # Rank by test statistic
  significant = "pval",
  p_cutoff = 0.05
)

# Display summary
cat("  LSEA summary (all characteristics):\n")
print(summary(lsea_all_result))

# Save the barplot
png(file.path(base_dir, "03.Enrichment/LSEA/LSEA_all_characteristics_barplot.png"),
  width = 1200, height = 1000, res = 150
)
print(lsea_all_result$static_barPlot)
dev.off()

cat("  ✓ LSEA (all characteristics) complete\n\n")


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

# Save the barplot
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

# Save the barplot
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

# ORA identifies enriched lipid classes among significant lipid species
# It classifies results into up-regulated, down-regulated, and non-significant

# ORA with all characteristics
cat("Running ORA for all characteristics...\n")
ora_all_result <- enrichment_ora(
  deSp_se,
  char = NULL, # NULL = all characteristics
  significant = "pval",
  p_cutoff = 0.05
)

# Display summary
cat("  ORA summary (all characteristics):\n")
print(summary(ora_all_result))

# Save the barplot (top 10 up- and down-regulated terms)
png(file.path(base_dir, "03.Enrichment/ORA/ORA_all_characteristics_barplot.png"),
  width = 1200, height = 1000, res = 150
)
print(ora_all_result$static_barPlot)
dev.off()

cat("  ✓ ORA (all characteristics) complete\n\n")


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

# Save the barplot
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

# Save the barplot
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
