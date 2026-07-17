##########################################################################################
# OXYLIPID CROSS-SECTIONAL ANALYSIS: MILD vs SEVERE, AT EACH OF D0/D3/D10/D60
##########################################################################################
#
# GOAL: for each timepoint independently, test whether individual oxylipins (and their
# parent-PUFA families/classes) differ between mild and severe dengue.
#
# DATA SOURCE: data/lipidsig_datasets/Oxylipines_raw/Oxylipid_lipid_abundance_data_*.tsv +
# group_information_table_oxylipid_*.tsv, built by code/Oxylipid_lipidomics_pretreatment.R
# from the raw (pre-ratio) abundances in data/raw_data/oxylipines.xlsx. There is no healthy
# arm for oxylipids (healthy donors are single-draw, no D60 sample), so this is mild vs
# severe only -- the oxylipid analogue of 02/03/04_posthoc_*.R, looped over 4 timepoints
# instead of 3 group pairs.
#
# WHY NOT LipidSigR's as_summarized_experiment() / deSp_twoGroup()?
# Oxylipins (HETEs, HODEs, prostaglandins, ...) are fatty-acid oxidation products, not
# glycerophospholipids/sphingolipids, so rgoslin::parseLipidNames() can't annotate them.
# A hand-built substitute annotation table was tried and does not work: LipidSigR's
# as_summarized_experiment() unconditionally runs its internal lipid_annotation(), which
# performs real LION/ChEBI/LIPID MAPS ontology lookups and nomenclature parsing that
# crash on non-rgoslin names (confirmed empirically -- a row-count mismatch when building
# the SummarizedExperiment). So this script reimplements the same statistical design
# directly in plain R: percentage normalization + log10 transform (matching
# data_process(normalization='Percentage', transform='log10')), t-test/Wilcoxon DE with
# BH-FDR (matching deSp_twoGroup(test='t-test')), and a Fisher's-exact class/category
# enrichment check (the same pattern 05_longitudinal_analysis.R already uses in place of
# LipidSigR's enrichment_ora() when no deSp_se object is available).
##########################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggsignif)
  library(ComplexHeatmap)
  library(circlize)
})

base_dir <- file.path(getwd(), "analysis", "Oxylipids", "PostHoc_Mild_vs_Severe")
data_dir <- "data/lipidsig_datasets/Oxylipines_raw"

class_annotation <- read_tsv(file.path(data_dir, "oxylipin_class_annotation.tsv"), show_col_types = FALSE)

timepoints <- c("D0", "D03", "D10", "D60")
names(timepoints) <- c("D0", "D3", "D10", "D60")  # display label -> file suffix

REF_GROUP <- "mild"  # log2FC / test direction: severe relative to mild

##########################################################################################
# Helper functions (shared across timepoints)
##########################################################################################

# Drop features detected (non-zero) in fewer than exclude_missing_pct% of samples --
# matches data_process(exclude_missing=TRUE, exclude_missing_pct=70) elsewhere in this
# project. Needed for the expanded 47-compound panel: several resolvins/protectins/minor
# prostaglandins are undetectable in most samples at a given timepoint, which would
# otherwise leave a constant (all-floor-value) row after normalize_log10's zero-replacement
# and break prcomp()'s per-feature scaling.
filter_low_prevalence <- function(abund_mat, exclude_missing_pct = 70) {
  detected_pct <- 100 * rowMeans(abund_mat > 0 & !is.na(abund_mat))
  keep <- !is.na(detected_pct) & detected_pct >= exclude_missing_pct
  if (any(!keep)) {
    cat("  Dropping", sum(!keep), "feature(s) detected in <", exclude_missing_pct, "% of samples:",
        paste(rownames(abund_mat)[!keep], collapse = ", "), "\n")
  }
  abund_mat[keep, , drop = FALSE]
}

# Percentage normalization (of total oxylipin signal per sample) + log10, with
# minimum-value replacement for zeros -- matches data_process(normalization='Percentage',
# replace_na_method='min', replace_na_method_ref=0.5, transform='log10') semantics.
normalize_log10 <- function(abund_mat) {
  pct <- sweep(abund_mat, 2, colSums(abund_mat, na.rm = TRUE), FUN = "/") * 100
  pct_filled <- t(apply(pct, 1, function(row) {
    nonzero_min <- suppressWarnings(min(row[row > 0], na.rm = TRUE))
    if (!is.finite(nonzero_min)) nonzero_min <- 1e-6
    row[is.na(row) | row <= 0] <- nonzero_min * 0.5
    row
  }))
  dimnames(pct_filled) <- dimnames(pct)
  log10(pct_filled)
}

fisher_enrichment <- function(sig_features, background_features, char_df, char_col, label) {
  if (length(sig_features) == 0) {
    return(data.frame(level = character(0), group = character(0), n_in_set = integer(0),
                       n_total_in_group = integer(0), odds_ratio = numeric(0),
                       p_value = numeric(0), fdr = numeric(0)))
  }
  bg <- char_df %>% filter(feature %in% background_features)
  bg %>%
    distinct(.data[[char_col]]) %>%
    pull(1) %>%
    lapply(function(grp) {
      in_group <- bg[[char_col]] == grp
      in_set <- bg$feature %in% sig_features
      tbl <- table(factor(in_set, levels = c(FALSE, TRUE)), factor(in_group, levels = c(FALSE, TRUE)))
      ft <- tryCatch(fisher.test(tbl), error = function(e) list(estimate = NA_real_, p.value = NA_real_))
      data.frame(
        level = label, group = grp,
        n_in_set = sum(in_set & in_group), n_total_in_group = sum(in_group),
        odds_ratio = unname(ft$estimate), p_value = ft$p.value
      )
    }) %>%
    bind_rows() %>%
    mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    arrange(fdr)
}

##########################################################################################
# Per-timepoint analysis
##########################################################################################

run_timepoint <- function(tp_label, tp_suffix) {
  cat("\n========================================\n")
  cat("TIMEPOINT:", tp_label, "(mild vs severe)\n")
  cat("========================================\n")

  out_dir <- file.path(base_dir, tp_label)
  dir.create(file.path(out_dir, "01.Profiling"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "02.DiffExp/Visualizations"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "02.DiffExp/Individual_lipid_boxplots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "03.Enrichment"), recursive = TRUE, showWarnings = FALSE)

  abundance <- read_tsv(file.path(data_dir, sprintf("Oxylipid_lipid_abundance_data_%s.tsv", tp_suffix)), show_col_types = FALSE)
  group_info <- read_tsv(file.path(data_dir, sprintf("group_information_table_oxylipid_%s.tsv", tp_suffix)), show_col_types = FALSE)
  group_info <- group_info %>% mutate(group = factor(group, levels = c(REF_GROUP, setdiff(unique(group), REF_GROUP))))

  cat("  ", nrow(abundance), "features x", nrow(group_info), "patients (",
      sum(group_info$group == "mild"), "mild,", sum(group_info$group == "severe"), "severe)\n")

  abund_mat <- as.matrix(abundance[, group_info$sample_name])
  rownames(abund_mat) <- abundance$feature
  abund_mat <- filter_low_prevalence(abund_mat)

  # -------------------- 1. Profiling --------------------
  processed <- normalize_log10(abund_mat)

  qc_df <- as.data.frame(processed) %>%
    tibble::rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample_name", values_to = "value") %>%
    left_join(group_info %>% select(sample_name, group), by = "sample_name")
  p_qc <- ggplot(qc_df, aes(x = sample_name, y = value, fill = group)) +
    geom_boxplot(outlier.size = 0.5) +
    labs(title = sprintf("%s: log10(%% of total oxylipin signal) per sample", tp_label), x = NULL, y = "log10(%)") +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
  png(file.path(out_dir, "01.Profiling/BoxPlot_normalized_by_sample.png"), width = 1800, height = 900, res = 150)
  print(p_qc); dev.off()

  pca <- prcomp(t(processed), scale. = TRUE, center = TRUE)
  var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  pca_df <- as.data.frame(pca$x[, 1:2]) %>%
    tibble::rownames_to_column("sample_name") %>%
    left_join(group_info %>% select(sample_name, group), by = "sample_name")
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 2.5) +
    labs(title = sprintf("%s: PCA (oxylipins)", tp_label),
         x = sprintf("PC1 (%.1f%%)", var_explained[1]), y = sprintf("PC2 (%.1f%%)", var_explained[2])) +
    theme_bw()
  png(file.path(out_dir, "01.Profiling/PCA.png"), width = 1200, height = 1000, res = 150)
  print(p_pca); dev.off()
  write_tsv(pca_df, file.path(out_dir, "01.Profiling/table_pca_scores.tsv"))

  sample_cor <- cor(processed, method = "spearman")
  ha <- HeatmapAnnotation(group = group_info$group[match(colnames(sample_cor), group_info$sample_name)],
                           col = list(group = c(mild = "#4DAF4A", severe = "#E41A1C")))
  png(file.path(out_dir, "01.Profiling/Heatmap_sample_correlation.png"), width = 1400, height = 1200, res = 150)
  draw(Heatmap(sample_cor, name = "Spearman r", top_annotation = ha,
               col = colorRamp2(c(min(sample_cor), 1), c("white", "#08519C")),
               column_title = sprintf("%s: sample-sample correlation", tp_label)))
  dev.off()

  # -------------------- 2. Differential expression --------------------
  mild_samples <- group_info$sample_name[group_info$group == "mild"]
  severe_samples <- group_info$sample_name[group_info$group == "severe"]

  de_results <- lapply(rownames(processed), function(f) {
    v_mild <- processed[f, mild_samples]
    v_severe <- processed[f, severe_samples]
    t_res <- tryCatch(t.test(v_severe, v_mild), error = function(e) NULL)
    w_res <- tryCatch(wilcox.test(v_severe, v_mild), error = function(e) NULL)
    log2fc <- (mean(v_severe, na.rm = TRUE) - mean(v_mild, na.rm = TRUE)) / log10(2)
    data.frame(
      feature = f, log2FC_severe_vs_mild = log2fc,
      p_ttest = if (is.null(t_res)) NA_real_ else t_res$p.value,
      p_wilcoxon = if (is.null(w_res)) NA_real_ else w_res$p.value
    )
  }) %>% bind_rows() %>%
    mutate(
      padj_ttest = p.adjust(p_ttest, method = "BH"),
      padj_wilcoxon = p.adjust(p_wilcoxon, method = "BH")
    ) %>%
    arrange(padj_ttest)

  write_tsv(de_results, file.path(out_dir, "02.DiffExp/DE_results_mild_vs_severe.tsv"))
  sig_results <- de_results %>% filter(padj_ttest < 0.05)
  write_tsv(sig_results, file.path(out_dir, "02.DiffExp/DE_significant_mild_vs_severe.tsv"))
  cat("  Significant (t-test padj<0.05):", nrow(sig_results), "/", nrow(de_results), "\n")

  # Volcano plot
  volcano_df <- de_results %>% mutate(sig = padj_ttest < 0.05)
  p_volcano <- ggplot(volcano_df, aes(x = log2FC_severe_vs_mild, y = -log10(p_ttest), color = sig)) +
    geom_point(size = 2) +
    ggrepel::geom_text_repel(data = filter(volcano_df, sig), aes(label = feature), size = 3, color = "black", max.overlaps = 20) +
    scale_color_manual(values = c(`TRUE` = "#E41A1C", `FALSE` = "grey60")) +
    labs(title = sprintf("%s: Volcano (severe vs mild)", tp_label), x = "log2FC (severe vs mild)", y = "-log10(p, t-test)") +
    theme_bw() + theme(legend.position = "none")
  png(file.path(out_dir, "02.DiffExp/Visualizations/Volcano.png"), width = 1400, height = 1200, res = 150)
  print(p_volcano); dev.off()

  # Boxplots for all 24 lipids (small enough to show every one, not just top hits)
  box_df <- as.data.frame(processed) %>%
    tibble::rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample_name", values_to = "value") %>%
    left_join(group_info %>% select(sample_name, group), by = "sample_name") %>%
    left_join(de_results %>% select(feature, padj_ttest), by = "feature") %>%
    mutate(sig_label = case_when(
      padj_ttest < 0.001 ~ "***", padj_ttest < 0.01 ~ "**", padj_ttest < 0.05 ~ "*", TRUE ~ "ns"
    ))

  feat_order <- de_results$feature
  p_box <- ggplot(box_df %>% mutate(feature = factor(feature, levels = feat_order)),
                   aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.5) +
    facet_wrap(~feature, scales = "free_y", ncol = 6) +
    scale_fill_manual(values = c(mild = "#4DAF4A", severe = "#E41A1C")) +
    labs(title = sprintf("%s: all oxylipins, mild vs severe (sorted by significance)", tp_label), x = NULL, y = "log10(%)") +
    theme_bw() + theme(legend.position = "none", strip.text = element_text(size = 7))
  png(file.path(out_dir, "02.DiffExp/Individual_lipid_boxplots/All_lipids.png"), width = 2400, height = 1800, res = 150)
  print(p_box); dev.off()

  # -------------------- 3. Class/category enrichment (Fisher's exact) --------------------
  all_features <- rownames(processed)
  enrichment_results <- bind_rows(
    fisher_enrichment(sig_results$feature, all_features, class_annotation, "class", "class"),
    fisher_enrichment(sig_results$feature, all_features, class_annotation, "family", "family")
  )
  write_tsv(enrichment_results, file.path(out_dir, "03.Enrichment/Class_family_enrichment.tsv"))
  cat("  Enrichment: tested", n_distinct(enrichment_results$group), "class/family group(s) (descriptive -- 24 features gives limited power)\n")

  invisible(list(timepoint = tp_label, n_patients = nrow(group_info), n_significant = nrow(sig_results)))
}

results_summary <- lapply(names(timepoints), function(tp) run_timepoint(tp, timepoints[[tp]]))

cat("\n========================================\n")
cat("CROSS-SECTIONAL OXYLIPID ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("Output directory:", base_dir, "\n")
for (r in results_summary) {
  cat(sprintf("  %s: %d patients, %d significant lipids (mild vs severe, t-test padj<0.05)\n", r$timepoint, r$n_patients, r$n_significant))
}
