##########################################################################################
# PLASMA LONGITUDINAL LMM: NORMALIZATION SENSITIVITY CHECK (PERCENTAGE vs PQN)
##########################################################################################
#
# Quick DE-only check (no trajectory plots, PCA, class/category aggregation, or
# enrichment) mirroring code/05_longitudinal_analysis.R's Sections 1-4: does the
# longitudinal LMM result (D0->D3->D10->D60, mild vs severe, the plasma pipeline's
# headline finding of 23/119 interaction-significant lipids) change under PQN
# normalization instead of the project-standard Percentage?
##########################################################################################

suppressPackageStartupMessages({
  library(LipidSigR)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(rgoslin)
  library(lme4)
  library(stringr)
})

setwd("E:/Dengue_lipidomics")

base_dir <- "E:/Dengue_lipidomics/analysis/Normalization_sensitivity/Longitudinal"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("PLASMA LONGITUDINAL LMM: NORMALIZATION SENSITIVITY CHECK\n")
cat("========================================\n")

cat("\n=== 1. DATA LOADING ===\n")

abund_files <- list(
  D0  = "data/lipidsig_datasets/Longitudinal_matched/Lipid_abundance_data_D0.tsv",
  D3  = "data/lipidsig_datasets/Longitudinal_matched/Lipid_abundance_data_D03.tsv",
  D10 = "data/lipidsig_datasets/Longitudinal_matched/Lipid_abundance_data_D10.tsv",
  D60 = "data/lipidsig_datasets/Longitudinal_matched/Lipid_abundance_data_D60.tsv"
)
group_files <- list(
  D0  = "data/lipidsig_datasets/Longitudinal_matched/group_information_table_D0.xlsx",
  D3  = "data/lipidsig_datasets/Longitudinal_matched/group_information_table_D03.xlsx",
  D10 = "data/lipidsig_datasets/Longitudinal_matched/group_information_table_D10.xlsx",
  D60 = "data/lipidsig_datasets/Longitudinal_matched/group_information_table_D60.xlsx"
)

abund_list <- lapply(abund_files, read_tsv, show_col_types = FALSE)
stopifnot(all(vapply(abund_list, function(x) identical(x$feature, abund_list$D0$feature), logical(1))))

abund_combined <- abund_list$D0
for (nm in c("D3", "D10", "D60")) {
  abund_combined <- dplyr::full_join(abund_combined, abund_list[[nm]], by = "feature")
}

group_combined <- bind_rows(lapply(group_files, read_excel))

mat <- as.matrix(abund_combined[, -1])
rownames(mat) <- abund_combined$feature
feat_const <- apply(mat, 1, function(x) length(unique(x)) <= 1)
samp_const <- apply(mat, 2, function(x) length(unique(x)) <= 1)

abund_clean <- abund_combined %>% dplyr::filter(!feat_const)
keep_samples <- colnames(mat)[!samp_const]
abund_clean <- abund_clean[, c("feature", keep_samples)]
group_clean <- group_combined %>% dplyr::filter(sample_name %in% keep_samples)
# Same positional-order requirement as the D0 three-group check
group_clean <- group_clean[match(colnames(abund_clean)[-1], group_clean$sample_name), ]
stopifnot(identical(group_clean$sample_name, colnames(abund_clean)[-1]))

cat("=== 2. LIPID ANNOTATION ===\n")

parse_lipid <- rgoslin::parseLipidNames(lipidNames = abund_clean$feature)
recognized_lipids <- parse_lipid$Original.Name[which(parse_lipid$Grammar != "NOT_PARSEABLE")]
abundance_filtered <- abund_clean %>% dplyr::filter(feature %in% recognized_lipids)
goslin_annotation <- parse_lipid %>% dplyr::filter(Original.Name %in% recognized_lipids)
cat("  Recognized lipids:", length(recognized_lipids), "/", nrow(abund_clean), "\n")

se <- as_summarized_experiment(
  abundance_filtered, goslin_annotation,
  group_info = group_clean, se_type = "profiling", paired_sample = NULL
)

sample_meta <- group_clean %>%
  mutate(
    timepoint = case_when(
      str_detect(sample_name, "D60$") ~ "D60",
      str_detect(sample_name, "D10$") ~ "D10",
      str_detect(sample_name, "D03$") ~ "D3",
      str_detect(sample_name, "D0$")  ~ "D0",
      TRUE ~ NA_character_
    ),
    patient_id = str_remove(sample_name, "D(60|10|03|0)$")
  ) %>%
  rename(severity = group) %>%
  select(sample_name, patient_id, timepoint, severity)

day_of_fever <- read_tsv("data/sick_patients_day_of_fever.tsv", show_col_types = FALSE)
sample_meta <- sample_meta %>% left_join(day_of_fever, by = c("patient_id" = "patient"))

fit_one_feature <- function(f, df) {
  d <- df %>% filter(feature == f)
  tryCatch({
    suppressMessages(withCallingHandlers({
      m_full   <- lmer(value ~ timepoint * severity + day_of_fever + (1 | patient_id), data = d, REML = FALSE)
      m_noint  <- lmer(value ~ timepoint + severity   + day_of_fever + (1 | patient_id), data = d, REML = FALSE)
      m_notime <- lmer(value ~ severity                + day_of_fever + (1 | patient_id), data = d, REML = FALSE)
      m_nosev  <- lmer(value ~ timepoint                + day_of_fever + (1 | patient_id), data = d, REML = FALSE)

      p_interaction <- anova(m_full, m_noint)[["Pr(>Chisq)"]][2]
      p_time        <- anova(m_noint, m_notime)[["Pr(>Chisq)"]][2]
      p_severity    <- anova(m_noint, m_nosev)[["Pr(>Chisq)"]][2]

      data.frame(feature = f, p_time = p_time, p_severity = p_severity, p_interaction = p_interaction, error = NA_character_)
    }, warning = function(w) invokeRestart("muffleWarning")))
  }, error = function(e) {
    data.frame(feature = f, p_time = NA_real_, p_severity = NA_real_, p_interaction = NA_real_, error = conditionMessage(e))
  })
}

run_lmm_check <- function(normalization_method) {
  cat("\n---", normalization_method, "---\n")
  processed_se <- tryCatch(
    data_process(
      se,
      exclude_missing = TRUE, exclude_missing_pct = 70,
      replace_na_method = "min", replace_na_method_ref = 0.5,
      normalization = normalization_method, transform = "log10"
    ),
    error = function(e) {
      cat("  SKIPPED (data_process) -", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(processed_se)) return(NULL)

  processed_abund <- processed_se@metadata$processed_abund
  cat("  Processed lipids:", nrow(processed_abund), "\n")

  long_df <- processed_abund %>%
    pivot_longer(-feature, names_to = "sample_name", values_to = "value") %>%
    inner_join(sample_meta, by = "sample_name") %>%
    mutate(timepoint = factor(timepoint, levels = c("D0", "D3", "D10", "D60")))

  features <- unique(long_df$feature)
  cat("  Fitting LMMs for", length(features), "lipids...\n")
  res <- bind_rows(lapply(features, fit_one_feature, df = long_df))

  n_failed <- sum(!is.na(res$error))
  if (n_failed > 0) cat("  ", n_failed, "lipid(s) failed to fit\n")

  res <- res %>%
    mutate(
      fdr_time = p.adjust(p_time, method = "BH"),
      fdr_interaction = p.adjust(p_interaction, method = "BH")
    ) %>%
    dplyr::select(feature, p_time, fdr_time, p_interaction, fdr_interaction) %>%
    dplyr::rename_with(~ paste0(., "_", normalization_method), -feature)

  cat("  Significant TIME effect (FDR<0.05):", sum(res[[paste0("fdr_time_", normalization_method)]] < 0.05, na.rm = TRUE), "/", nrow(res), "\n")
  cat("  Significant INTERACTION (FDR<0.05):", sum(res[[paste0("fdr_interaction_", normalization_method)]] < 0.05, na.rm = TRUE), "/", nrow(res), "\n")

  res
}

cat("\n=== 3. RUNNING LMM UNDER EACH NORMALIZATION METHOD ===\n")
results <- lapply(c("Percentage", "PQN"), run_lmm_check)
names(results) <- c("Percentage", "PQN")
results <- results[!sapply(results, is.null)]

comparison <- Reduce(function(x, y) dplyr::full_join(x, y, by = "feature"), results)
write_tsv(comparison, file.path(base_dir, "Longitudinal_normalization_comparison.tsv"))

cat("\n=== 4. INTERACTION-SIGNIFICANT LIPID OVERLAP ===\n")
sig_int_pct <- comparison$feature[comparison$fdr_interaction_Percentage < 0.05 & !is.na(comparison$fdr_interaction_Percentage)]
sig_int_pqn <- comparison$feature[comparison$fdr_interaction_PQN < 0.05 & !is.na(comparison$fdr_interaction_PQN)]
cat("  Percentage interaction-significant:", length(sig_int_pct), "\n")
cat("  PQN interaction-significant:", length(sig_int_pqn), "\n")
cat("  In both:", length(intersect(sig_int_pct, sig_int_pqn)), "->", paste(intersect(sig_int_pct, sig_int_pqn), collapse = ", "), "\n")
cat("  Percentage only:", length(setdiff(sig_int_pct, sig_int_pqn)), "->", paste(setdiff(sig_int_pct, sig_int_pqn), collapse = ", "), "\n")
cat("  PQN only:", length(setdiff(sig_int_pqn, sig_int_pct)), "->", paste(setdiff(sig_int_pqn, sig_int_pct), collapse = ", "), "\n")

cat("\n=== TIME-EFFECT LIPID OVERLAP ===\n")
sig_time_pct <- comparison$feature[comparison$fdr_time_Percentage < 0.05 & !is.na(comparison$fdr_time_Percentage)]
sig_time_pqn <- comparison$feature[comparison$fdr_time_PQN < 0.05 & !is.na(comparison$fdr_time_PQN)]
cat("  Percentage time-significant:", length(sig_time_pct), "\n")
cat("  PQN time-significant:", length(sig_time_pqn), "\n")
cat("  In both:", length(intersect(sig_time_pct, sig_time_pqn)), "\n")

cat("\nOutputs written to:", base_dir, "\n")
