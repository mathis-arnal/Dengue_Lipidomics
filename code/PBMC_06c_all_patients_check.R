##########################################################################################
# PBMC HEALTHY vs DISEASED: DOES THE RESULT HOLD WITH ALL PATIENTS (NOT JUST THE
# 40-PATIENT PLASMA-MATCHED COHORT)?
##########################################################################################
#
# PBMC_02/03/04/06 restrict to the canonical 40-patient cohort (20 healthy + 20 sick)
# because it's the one matched 1:1 against the plasma dataset's patient list. But PBMC
# itself has abundance data for more patients than that -- group_information_table_
# healthy_vs_sick_patients_D0_all.tsv covers 63 patients (20 healthy + 22 mild + 21
# severe), built by PBMC_01_analysis_three_groups_all.R. This script re-derives the
# healthy-vs-diseased pooling from that larger cohort and re-runs the same DE test to
# check whether the 34-patient (matched-to-plasma) result generalizes to the larger,
# PBMC-only sample.
##########################################################################################

suppressPackageStartupMessages({
  library(LipidSigR)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(rgoslin)
})

source("code/lipidomics_functions.R")

setwd("E:/Dengue_lipidomics")

base_dir <- "E:/Dengue_lipidomics/analysis/PBMC/PostHoc_Healthy_vs_Diseased/D0/All_patients_check"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("HEALTHY vs DISEASED: ALL-PATIENTS CHECK\n")
cat("========================================\n")

cat("\nDeriving healthy-vs-diseased group information from the ALL (63-patient) cohort...\n")
all_group_info <- read_tsv("data/PBMC/group_information_table_healthy_vs_sick_patients_D0_all.tsv", show_col_types = FALSE)
diseased_all_group_info <- all_group_info %>%
  mutate(group = ifelse(group %in% c("mild", "severe"), "diseased", group))
write_tsv(diseased_all_group_info, "data/PBMC/group_information_table_healthy_vs_diseased_D0_all.tsv")
cat("  Healthy:", sum(diseased_all_group_info$group == "healthy"),
    " Diseased:", sum(diseased_all_group_info$group == "diseased"), "\n")

loaded <- load_and_match_samples(
  "data/PBMC/group_information_table_healthy_vs_diseased_D0_all.tsv",
  "data/PBMC/PBMC_Lipid_abundance_data_D0.tsv"
)
abundance <- loaded$abundance
group_info <- loaded$group_info

annotated <- clean_and_annotate_lipids(abundance)

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
      df <- as.data.frame(res_list$all_deSp_result) %>%
        dplyr::select(feature, log2FC, pval, padj) %>%
        dplyr::rename_with(~ paste0(., "_", normalization_method), c(log2FC, pval, padj))
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

results <- lapply(c("Percentage", "PQN", "none"), run_de)
names(results) <- c("Percentage", "PQN", "none")
results <- results[!sapply(results, is.null)]

comparison <- Reduce(function(x, y) dplyr::full_join(x, y, by = "feature"), results)
comparison <- comparison %>% dplyr::arrange(pval_Percentage)
write_tsv(comparison, file.path(base_dir, "all_patients_comparison.tsv"))

cat("\nOutputs written to:", base_dir, "\n")
