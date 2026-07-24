##########################################################################################
# PBMC HEALTHY vs DISEASED: NORMALIZATION SENSITIVITY CHECK
##########################################################################################
#
# Does the choice of data_process() normalization method change the Healthy-vs-Diseased
# DE result? Percentage normalization (the project standard, used by PBMC_06) forces
# every sample's lipidome to sum to 100%, which risks distorting individual-lipid results
# specifically because PBMC_05_total_abundance_check.R already found that total lipid
# abundance itself differs significantly by severity (Kruskal-Wallis p=0.021) -- a real
# compositional-closure concern for biomarker discovery, not just a theoretical one.
#
# Compares three methods on the SAME data/samples/features:
#   - Percentage : the project standard (what PBMC_06 uses)
#   - PQN        : Probabilistic Quotient Normalization -- more robust to the "one real
#                  change forces every other feature to shift" problem than a fixed-sum
#                  method; already used in this project's original oxylipid prototype
#                  (LipidSigR/Profiling.R)
#   - none       : no normalization at all -- lets any total-abundance-driven signal
#                  through uncorrected, as a sanity-check bookend
# ("Sum" is skipped: statistically identical to Percentage after log-transform, since
# dividing by total either as a fraction or x100 is just a constant additive shift in
# log-space that cannot change any p-value. "Quantile"/"Median" skipped as lower-priority
# per the same reasoning discussed for this check.)
#
# Scope: Healthy vs Diseased only (PBMC_06's comparison) -- the best-powered PBMC
# comparison run so far (34 patients, 8 nominal hits under Percentage).
##########################################################################################

suppressPackageStartupMessages({
  library(LipidSigR)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(rgoslin)
  library(ggplot2)
})

source("code/lipidomics_functions.R")

setwd("E:/Dengue_lipidomics")

base_dir <- "E:/Dengue_lipidomics/analysis/PBMC/PostHoc_Healthy_vs_Diseased/D0/Normalization_sensitivity"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("NORMALIZATION SENSITIVITY CHECK: HEALTHY vs DISEASED\n")
cat("========================================\n")

cat("\nLoading group information and abundance data...\n")
loaded <- load_and_match_samples(
  "data/PBMC/group_information_table_healthy_vs_diseased_D0.tsv",
  "data/PBMC/PBMC_Lipid_abundance_data_D0.tsv"
)
abundance <- loaded$abundance
group_info <- loaded$group_info

annotated <- clean_and_annotate_lipids(abundance)

# PQN requires the 'Rcpm' package, which LipidSigR declares as a dependency but which
# has no resolved renv.lock entry and isn't on CRAN -- a genuine environment gap, not
# something to silently patch by installing an unverified package from an unclear
# source. Wrapped in tryCatch so Percentage/none still run and produce a result even if
# PQN is unavailable in this environment.
run_de <- function(normalization_method) {
  cat("\n---", normalization_method, "---\n")
  tryCatch(
    {
      se_result <- build_processed_se(
        annotated$abundance_filtered, annotated$goslin_annotation, group_info,
        se_type = "de_two", exclude_missing_pct = 70, paired_sample = FALSE,
        normalization = normalization_method
      )
      deSp_se <- deSp_twoGroup(
        se_result$processed_se,
        ref_group = "healthy", test = "t-test",
        significant = "pval", p_cutoff = 0.05, FC_cutoff = 1, transform = "log10"
      )
      res_list <- extract_summarized_experiment(deSp_se)
      # all_deSp_result (not sig_deSp_result) so every tested lipid is comparable across methods
      df <- as.data.frame(res_list$all_deSp_result) %>%
        dplyr::select(feature, log2FC, pval, padj) %>%
        dplyr::rename_with(~ paste0(., "_", normalization_method), c(log2FC, pval, padj))
      cat("  Lipids tested:", nrow(df), " | nominal p<0.05:", sum(res_list$all_deSp_result$pval < 0.05, na.rm = TRUE),
          " | FDR<0.05:", sum(res_list$all_deSp_result$padj < 0.05, na.rm = TRUE), "\n")
      df
    },
    error = function(e) {
      cat("  SKIPPED -", conditionMessage(e), "\n")
      NULL
    }
  )
}

results <- lapply(c("Percentage", "PQN", "none"), run_de)
names(results) <- c("Percentage", "PQN", "none")
unavailable <- names(results)[sapply(results, is.null)]
if (length(unavailable) > 0) {
  cat("\nSkipped (unavailable in this environment):", paste(unavailable, collapse = ", "), "\n")
}
results <- results[!sapply(results, is.null)]

# Merge into one wide comparison table (inner join -- same features tested regardless of
# normalization method, since normalization doesn't change which lipids pass exclude_missing_pct)
comparison <- Reduce(function(x, y) dplyr::full_join(x, y, by = "feature"), results)
comparison <- comparison %>% dplyr::arrange(pval_Percentage)
write_tsv(comparison, file.path(base_dir, "normalization_comparison.tsv"))

cat("\n=== Nominal significance (p<0.05) overlap across methods ===\n")
methods_available <- names(results)
sig_sets <- lapply(methods_available, function(m) {
  comparison$feature[comparison[[paste0("pval_", m)]] < 0.05 & !is.na(comparison[[paste0("pval_", m)]])]
})
names(sig_sets) <- methods_available
for (m in names(sig_sets)) cat(" ", m, ":", length(sig_sets[[m]]), "->", paste(sig_sets[[m]], collapse = ", "), "\n")

if (length(sig_sets) > 1) {
  cat("\n  In ALL", length(sig_sets), "methods:", paste(Reduce(intersect, sig_sets), collapse = ", "), "\n")
  for (m in names(sig_sets)) {
    others <- setdiff(names(sig_sets), m)
    only_m <- setdiff(sig_sets[[m]], Reduce(union, sig_sets[others]))
    cat(" ", m, "only:", paste(only_m, collapse = ", "), "\n")
  }
}

# Rank correlation between methods' p-values -- are they measuring roughly the same
# thing (just rescaled), or genuinely disagreeing on which lipids stand out?
cat("\n=== Spearman correlation between methods' p-values (log10) ===\n")
p_mat <- comparison %>% dplyr::select(starts_with("pval_")) %>% dplyr::mutate(across(everything(), log10))
print(round(cor(p_mat, method = "spearman", use = "pairwise.complete.obs"), 3))

# Volcano-style comparison plot: log2FC and -log10(p) side by side per method
plot_df <- comparison %>%
  tidyr::pivot_longer(
    cols = -feature,
    names_to = c(".value", "method"),
    names_pattern = "(log2FC|pval|padj)_(.*)"
  )
p <- ggplot(plot_df, aes(x = log2FC, y = -log10(pval), color = padj < 0.05)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  facet_wrap(~method) +
  scale_color_manual(values = c(`TRUE` = "#E41A1C", `FALSE` = "grey40")) +
  labs(
    title = "Healthy vs Diseased: DE result by normalization method",
    subtitle = "Dashed line = nominal p=0.05 (not FDR); color = FDR<0.05",
    x = "log2FC (diseased vs healthy)", y = "-log10(p, t-test)", color = "FDR<0.05"
  ) +
  theme_bw()
ggsave(file.path(base_dir, "normalization_comparison_volcano.png"), p, width = 12, height = 5, dpi = 150)

cat("\nOutputs written to:", base_dir, "\n")
