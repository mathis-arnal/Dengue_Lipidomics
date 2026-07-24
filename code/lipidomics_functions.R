##########################################################################################
# SHARED LIPIDOMICS ANALYSIS FUNCTIONS
##########################################################################################
#
# Extracted from PBMC_01_analysis_three_groups.R / PBMC_01_analysis_three_groups_all.R /
# PBMC_02_posthoc_healthy_vs_mild.R / PBMC_03_posthoc_healthy_vs_severe.R /
# PBMC_04_posthoc_mild_vs_severe.R, which duplicated the same data-loading,
# rgoslin-annotation, SummarizedExperiment-construction, and plot-saving boilerplate
# near-verbatim. Source this file at the top of each analysis script:
#
#   source("code/lipidomics_functions.R")
#
# Deliberately NOT extracted into shared functions (kept inline per script, since they
# differ in ways that matter for correctness, not just cosmetically):
#   - dr_pca()/dr_tsne()/dr_umap()/heatmap_correlation() calls themselves -- cluster_num,
#     and which SE object (se vs processed_se) is passed to cross_sample_variability(),
#     differ between the three-group script (cluster_num=3) and the posthoc scripts
#     (cluster_num=2) in ways that reflect real analysis choices, not accidents.
#   - The DE + boxplot block -- ANOVA (multi-group, boxPlot_feature_multiGroup(),
#     `adj.P.Val` column) vs t-test (two-group, manual ggplot2+ggsignif boxplots with
#     per-comparison box colors, `padj` column) are different enough that forcing them
#     into one function would obscure real logic differences rather than remove
#     duplication.
##########################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(rgoslin)
})


#' Load abundance + group-info TSVs and restrict both to their common samples.
#'
#' Mirrors the "load group_info, load abundance, intersect sample names, filter both,
#' print diagnostics" block repeated at the top of every PBMC_0*.R script.
#'
#' @param group_info_path Path to a `sample_name/label_name/group/pair` TSV.
#' @param abundance_path Path to a `feature` + one-column-per-sample abundance TSV.
#' @return list(abundance, group_info), both filtered to matched_samples only.
load_and_match_samples <- function(group_info_path, abundance_path) {
  group_info <- read_tsv(group_info_path, show_col_types = FALSE)
  abundance <- read_tsv(abundance_path, show_col_types = FALSE)

  available_samples <- colnames(abundance)
  matched_samples <- intersect(group_info$sample_name, available_samples)

  cat("Samples in group_info:", nrow(group_info), "\n")
  cat("Samples in abundance:", length(available_samples) - 1, "\n") # minus 'feature' col
  cat("Matched samples:", length(matched_samples), "\n")
  unmatched <- setdiff(group_info$sample_name, available_samples)
  if (length(unmatched) > 0) {
    cat("Unmatched from group_info:", paste(unmatched, collapse = ", "), "\n")
  }
  cat("\n")

  abundance <- abundance %>% dplyr::select(feature, dplyr::all_of(matched_samples))
  group_info <- group_info %>% dplyr::filter(sample_name %in% matched_samples)

  cat("  Features:", nrow(abundance), "\n")
  cat("  Samples:", ncol(abundance) - 1, "\n")
  cat("  Total patient count:", nrow(group_info), "\n")
  for (grp in sort(unique(group_info$group))) {
    cat("  ", grp, "patient count:", sum(group_info$group == grp), "\n")
  }

  list(abundance = abundance, group_info = group_info)
}


#' Clean lipid feature names and annotate them with rgoslin, keeping only
#' recognized (parseable) lipids.
#'
#' @param abundance A `feature` + sample-columns data frame.
#' @return list(abundance_filtered, goslin_annotation), both restricted to
#'   rgoslin-recognized lipids only.
clean_and_annotate_lipids <- function(abundance) {
  cat("Parsing lipid nomenclature with rgoslin...\n")

  # (1,2) DG 32:0 -> DG 32:0
  abundance$feature <- gsub("\\(\\d+,\\d+\\) DG", "DG", abundance$feature)
  # Cer 34:1,O2 -> Cer 34:1;O2
  abundance$feature <- gsub(",", ";", abundance$feature)

  parse_lipid <- rgoslin::parseLipidNames(lipidNames = abundance$feature)
  recognized_lipids <- parse_lipid$Original.Name[which(parse_lipid$Grammar != "NOT_PARSEABLE")]

  abundance_filtered <- abundance %>% dplyr::filter(feature %in% recognized_lipids)
  goslin_annotation <- parse_lipid %>% dplyr::filter(Original.Name %in% recognized_lipids)

  cat("  Recognized lipids:", length(recognized_lipids), "/", nrow(abundance), "\n")
  cat("  Recognition rate:", round(100 * length(recognized_lipids) / nrow(abundance), 1), "%\n")

  list(abundance_filtered = abundance_filtered, goslin_annotation = goslin_annotation)
}


#' Build and process a LipidSigR SummarizedExperiment with this project's standard
#' data_process() parameters (min-imputation x0.5, Percentage normalization, log10
#' transform -- see CLAUDE.md "Architecture").
#'
#' @param se_type "de_multiple" (three-group) or "de_two" (posthoc pairwise).
#' @param exclude_missing_pct Passed through to data_process(); differs between
#'   scripts (70 for Healthy vs Mild, 50 for the sparser Healthy vs Severe / Mild vs
#'   Severe comparisons) so it is a required argument, not a hardcoded default.
#' @param normalization One of data_process()'s normalization options
#'   ("none","Percentage","PQN","Quantile","Sum","Median"). Defaults to "Percentage"
#'   (this project's standard convention) so existing callers are unaffected; exposed
#'   as a parameter for normalization-sensitivity checks.
#' @return list(se, processed_se)
build_processed_se <- function(abundance_filtered, goslin_annotation, group_info,
                                se_type, exclude_missing_pct, paired_sample = FALSE,
                                normalization = "Percentage") {
  # as_summarized_experiment()'s internal .check_group_info() requires group_info to
  # have EXACTLY the columns matching se_type: sample_name/label_name/group/pair for
  # "de_two", but sample_name/label_name/group only (no pair) for "de_multiple" --
  # drop `pair` here rather than pushing this LipidSigR quirk onto every caller.
  if (se_type == "de_multiple" && "pair" %in% colnames(group_info)) {
    group_info <- group_info %>% dplyr::select(-pair)
  }

  cat("Creating SummarizedExperiment object...\n")
  se <- as_summarized_experiment(
    abundance_filtered, goslin_annotation,
    group_info = group_info, se_type = se_type, paired_sample = paired_sample
  )
  cat("  ✓ SummarizedExperiment created\n")

  cat("Processing data (filtering, normalization, transformation)...\n")
  processed_se <- data_process(
    se,
    exclude_missing = TRUE, exclude_missing_pct = exclude_missing_pct,
    replace_na_method = "min", replace_na_method_ref = 0.5,
    normalization = normalization, transform = "log10"
  )
  cat("  ✓ Data processing complete\n")

  list(se = se, processed_se = processed_se)
}


#' Build a clinical-group-colored PCA scatter on top of an already-computed dr_pca()
#' result, with unsupervised k-means cluster kept as point shape and (by default) a
#' 95% confidence ellipse per clinical group.
#'
#' dr_pca(..., clustering="group_info") cannot be used directly to get a plot colored
#' by the true clinical group: its group_info extraction depends on metadata added by
#' a prior deSp_multiGroup()/deSp_twoGroup() call, not on the se_type used to build the
#' SE, so it returns NULL there and dr_pca() errors ("clustering select 'group_info' in
#' processed_se must provide group_info"). Joining group_info back onto the
#' kmeans-clustered PCA scores is the working alternative used here.
#'
#' @param result_pca Output of dr_pca() (must have been run with clustering="kmeans").
#' @param group_info sample_name/label_name/group[/pair] table for the same samples.
#' @param ellipses Draw a 95% normal-confidence ellipse per clinical group (via
#'   ggplot2::stat_ellipse) -- matches this project's existing convention of describing
#'   PCA plots by group-ellipse shape/overlap. Silently skipped by ggplot2 for any group
#'   with too few points to fit a covariance ellipse.
#' @return list(plot = <ggplot object>, plot_interactive = <plotly object>)
plot_pca_by_severity <- function(result_pca, group_info, ellipses = TRUE) {
  pca_data <- result_pca$pca_rotated_data %>%
    dplyr::rename(cluster = group) %>%
    dplyr::mutate(cluster = factor(cluster)) %>%
    dplyr::left_join(group_info %>% dplyr::select(sample_name, group), by = "sample_name")

  g <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = group, shape = cluster, text = sample_name)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme_bw() +
    ggplot2::labs(color = "Clinical group", shape = "K-means cluster (unsupervised)")

  if (ellipses) {
    g <- g + ggplot2::stat_ellipse(
      ggplot2::aes(group = group, color = group),
      type = "norm", level = 0.95, linewidth = 0.6, show.legend = FALSE
    )
  }

  g_interactive <- plotly::ggplotly(g, tooltip = c("text", "colour", "shape"))

  list(plot = g, plot_interactive = g_interactive)
}


#' Create a set of nested output directories under base_dir.
#'
#' @param subdirs Character vector of relative subdirectory paths, e.g.
#'   c("01.Profiling/00.Data_quality", "02.DiffExp/01.ANOVA").
create_output_dirs <- function(base_dir, subdirs) {
  for (d in subdirs) {
    dir.create(file.path(base_dir, d), recursive = TRUE, showWarnings = FALSE)
  }
}


#' Save a static plot object to PNG. Silently does nothing if plot_obj is NULL
#' (mirrors the `if (!is.null(...))` guards used throughout these scripts).
save_static_png <- function(plot_obj, path, width = 1200, height = 1000, res = 150) {
  if (is.null(plot_obj)) {
    return(invisible(NULL))
  }
  png(path, width = width, height = height, res = res)
  print(plot_obj)
  dev.off()
}


#' Save an interactive (plotly/htmlwidget) plot object as .rds for later
#' embedding via readRDS() in an Rmd report. Silently does nothing if NULL
#' (LipidSigR returns NULL interactive plots in some edge cases, e.g. too few
#' groups/features for a given plot type).
save_interactive_rds <- function(interactive_obj, path) {
  if (is.null(interactive_obj)) {
    return(invisible(NULL))
  }
  saveRDS(interactive_obj, path)
}


#' Save an interactive plot as both a self-contained HTML widget and an .rds.
#' Used for the two-group DE Volcano/MA plots, which (unlike other interactive
#' plots in this pipeline) are also saved as standalone HTML.
save_interactive_html_and_rds <- function(interactive_obj, html_path, rds_path) {
  if (is.null(interactive_obj)) {
    return(invisible(NULL))
  }
  htmlwidgets::saveWidget(interactive_obj, html_path, selfcontained = TRUE)
  saveRDS(interactive_obj, rds_path)
}


#' Run enrichment_lsea() and enrichment_ora() over a set of lipid characteristics,
#' saving a barplot PNG for each and printing a summary.
#'
#' @param char_list Named list mapping a label (used in the log and filename) to
#'   the `char` value passed to enrichment_lsea()/enrichment_ora(), e.g.
#'   list(all_characteristics = NULL, by_class = "class", by_category = "Category").
#' @param width,height,res PNG dimensions -- 1200x1000@150 matches every LSEA/ORA
#'   plot in PBMC_01-04, so the default should not normally need overriding.
run_lsea_ora <- function(deSp_se, char_list, lsea_dir, ora_dir, p_cutoff = 0.05,
                          width = 1200, height = 1000, res = 150) {
  for (label in names(char_list)) {
    char_val <- char_list[[label]]

    cat("Running LSEA (", label, ")...\n")
    lsea_result <- enrichment_lsea(
      deSp_se,
      char = char_val, rank_by = "statistic",
      significant = "pval", p_cutoff = p_cutoff
    )
    cat("  LSEA summary (", label, "):\n")
    print(summary(lsea_result))
    save_static_png(
      lsea_result$static_barPlot,
      file.path(lsea_dir, paste0("LSEA_", label, "_barplot.png")),
      width = width, height = height, res = res
    )
    cat("  ✓ LSEA (", label, ") complete\n\n")

    cat("Running ORA (", label, ")...\n")
    ora_result <- enrichment_ora(
      deSp_se,
      char = char_val, significant = "pval", p_cutoff = p_cutoff
    )
    cat("  ORA summary (", label, "):\n")
    print(summary(ora_result))
    save_static_png(
      ora_result$static_barPlot,
      file.path(ora_dir, paste0("ORA_", label, "_barplot.png")),
      width = width, height = height, res = res
    )
    cat("  ✓ ORA (", label, ") complete\n\n")
  }
  invisible(NULL)
}
