##########################################################################################
# PLASMA D0 THREE-GROUP ANOVA: NORMALIZATION SENSITIVITY CHECK (PERCENTAGE vs PQN)
##########################################################################################
#
# Quick DE-only check (no profiling/plots/enrichment) mirroring the PBMC normalization-
# sensitivity script (PBMC_06b): does the D0 three-group ANOVA result (Healthy vs Mild
# vs Severe, the plasma pipeline's headline cross-sectional comparison) change under PQN
# normalization instead of the project-standard Percentage? Uses the same data source as
# code/01_analysis_three_groups.R.
##########################################################################################

suppressPackageStartupMessages({
  library(LipidSigR)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(rgoslin)
})

setwd("E:/Dengue_lipidomics")

base_dir <- "E:/Dengue_lipidomics/analysis/Normalization_sensitivity/D0_ThreeGroup"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("PLASMA D0 THREE-GROUP: NORMALIZATION SENSITIVITY CHECK\n")
cat("========================================\n")

cat("\nLoading group information and abundance data...\n")
group_info <- read_tsv("data/lipidsig_datasets/healthy_vs_sick_patients/group_information_table_healthy_vs_sick_patients_D0.tsv", show_col_types = FALSE)
abundance <- read_tsv("data/lipidsig_datasets/healthy_vs_sick_patients/healthy_sick_lipidomics.tsv", show_col_types = FALSE)
cat("  Groups:", paste(unique(group_info$group), collapse = ", "), " | Samples:", nrow(group_info), "\n")

cat("Parsing lipid nomenclature with rgoslin...\n")
parse_lipid <- rgoslin::parseLipidNames(lipidNames = abundance$feature)
recognized_lipids <- parse_lipid$Original.Name[which(parse_lipid$Grammar != "NOT_PARSEABLE")]
abundance_filtered <- abundance %>% dplyr::filter(feature %in% recognized_lipids)
goslin_annotation <- parse_lipid %>% dplyr::filter(Original.Name %in% recognized_lipids)
cat("  Recognized lipids:", length(recognized_lipids), "/", nrow(abundance), "\n")

# group_info and abundance have the same 40 samples but not in the same column order --
# as_summarized_experiment() requires an exact positional match, not just the same set.
group_info <- group_info[match(colnames(abundance_filtered)[-1], group_info$sample_name), ]
stopifnot(identical(group_info$sample_name, colnames(abundance_filtered)[-1]))

se <- as_summarized_experiment(
  abundance_filtered, goslin_annotation,
  group_info = group_info, se_type = "de_multiple", paired_sample = NULL
)

run_anova <- function(normalization_method) {
  cat("\n---", normalization_method, "---\n")
  tryCatch(
    {
      processed_se <- data_process(
        se,
        exclude_missing = TRUE, exclude_missing_pct = 70,
        replace_na_method = "min", replace_na_method_ref = 0.5,
        normalization = normalization_method, transform = "log10"
      )
      deSp_se <- deSp_multiGroup(
        processed_se, test = "One-way ANOVA", ref_group = "healthy",
        significant = "pval", p_cutoff = 0.05, transform = "log10"
      )
      res_list <- extract_summarized_experiment(deSp_se)
      df <- as.data.frame(res_list$all_deSp_result) %>%
        dplyr::select(feature, pval, padj) %>%
        dplyr::rename_with(~ paste0(., "_", normalization_method), c(pval, padj))
      cat("  Lipids tested:", nrow(df), " | nominal p<0.05:", sum(res_list$all_deSp_result$pval < 0.05, na.rm = TRUE),
          " | FDR<0.05:", sum(res_list$all_deSp_result$padj < 0.05, na.rm = TRUE), "\n")
      sig <- res_list$all_deSp_result$feature[res_list$all_deSp_result$pval < 0.05 & !is.na(res_list$all_deSp_result$pval)]
      cat("  Nominal hits:", paste(sig, collapse = ", "), "\n")
      df
    },
    error = function(e) {
      cat("  SKIPPED -", conditionMessage(e), "\n")
      NULL
    }
  )
}

results <- lapply(c("Percentage", "PQN"), run_anova)
names(results) <- c("Percentage", "PQN")
results <- results[!sapply(results, is.null)]

comparison <- Reduce(function(x, y) dplyr::full_join(x, y, by = "feature"), results)
comparison <- comparison %>% dplyr::arrange(pval_Percentage)
write_tsv(comparison, file.path(base_dir, "D0_threegroup_normalization_comparison.tsv"))

cat("\n=== Nominal significance (p<0.05) overlap ===\n")
sig_pct <- comparison$feature[comparison$pval_Percentage < 0.05 & !is.na(comparison$pval_Percentage)]
sig_pqn <- comparison$feature[comparison$pval_PQN < 0.05 & !is.na(comparison$pval_PQN)]
cat("  Percentage:", length(sig_pct), "->", paste(sig_pct, collapse = ", "), "\n")
cat("  PQN:", length(sig_pqn), "->", paste(sig_pqn, collapse = ", "), "\n")
cat("  In both:", paste(intersect(sig_pct, sig_pqn), collapse = ", "), "\n")
cat("  Percentage only:", paste(setdiff(sig_pct, sig_pqn), collapse = ", "), "\n")
cat("  PQN only:", paste(setdiff(sig_pqn, sig_pct), collapse = ", "), "\n")

cat("\nOutputs written to:", base_dir, "\n")
