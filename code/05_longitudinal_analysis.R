##########################################################################################
# LONGITUDINAL LIPIDOMICS ANALYSIS: D0 -> D3 -> D10 -> D60 (MILD vs SEVERE, PAIRED)
##########################################################################################
#
# GOAL: track how each patient's lipidome evolves across all 4 timepoints, rather than
# comparing groups cross-sectionally at a single timepoint.
#
# DATA SOURCE: data/lipidsig_datasets/Longitudinal_matched/Lipid_abundance_data_{D0,D03,D10,D60}.tsv
# + matching group_information_table_*.xlsx. This is the
# only dataset in the project where the SAME 43 patients (mild/severe, no healthy) are
# tracked with raw abundances at all 4 timepoints, in identical sample order/columns.
# There is no longitudinal healthy-control data anywhere in the project, so healthy
# patients are not part of this analysis.
#
# WHY NOT LipidSigR's deSp_multiGroup (ANOVA)?
# It only supports One-way ANOVA / Kruskal-Wallis, both of which assume independent
# samples. Feeding it D0/D3/D10/D60-as-groups would treat a patient's 4 repeated
# measurements as 4 independent observations (pseudoreplication), inflating false
# positives and ignoring each patient's own baseline. Instead, per-lipid statistics here
# use a linear mixed model with patient as a random intercept, which is the standard way
# to test change over time within subjects while still comparing mild vs severe.
#
# STATISTICAL MODEL (per lipid):
#   value ~ timepoint * severity + day_of_fever + (1 | patient_id)
# Timepoint is a 4-level FACTOR (not numeric days) because D0/D3/D10/D60 spacing is
# highly uneven (0, 3, 10, 60 days) -- treating it as linear would assume a trend that
# isn't justified by the design. day_of_fever (days of fever before hospitalization/
# sampling, from data/sick_patients_day_of_fever.tsv) is included as a covariate in every
# model below because D0 is enrollment day, not a fixed point in each patient's illness
# course -- two patients both sampled at "D0" can already be several days apart in actual
# disease progression, which would otherwise confound the severity/time comparisons.
# Models are fit by ML (REML=FALSE) so nested models can be compared via likelihood-ratio
# tests (LRT):
#   - p_time          : full-minus-interaction  vs  severity-only      -> does abundance
#                        change over time at all (averaged across groups)?
#   - p_severity      : full-minus-interaction  vs  timepoint-only     -> do groups differ
#                        on average (across time)?
#   - p_interaction   : full                     vs  full-minus-interaction -> does the
#                        trajectory shape differ between mild and severe?
# P-values for each effect are BH-FDR-adjusted across all tested lipids.
#
# LipidSigR is still used for what it's good at here: rgoslin lipid annotation, QC,
# normalization (Percentage + log10), and PCA, which are all per-sample operations that
# don't assume group independence. Its group-comparison DE/enrichment functions
# (deSp_multiGroup, LSEA, ORA) are NOT used, for the reason above.
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
  library(emmeans)
  library(stringr)
  library(ggplot2)
})

setwd("E:/PBMC_lipidomics")

base_dir <- file.path(getwd(), "analysis", "Longitudinal")
dir.create(file.path(base_dir, "00.Data_quality"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "01.Trajectories/PCA"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "02.Mixed_models"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "02.Mixed_models/Trajectory_plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "03.Aggregated/Class/Trajectory_plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "03.Aggregated/Category/Trajectory_plots"), recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("LONGITUDINAL LIPIDOMICS ANALYSIS (D0/D3/D10/D60)\n")
cat("Groups: Mild, Severe (paired, no healthy longitudinal data)\n")
cat("========================================\n")


##########################################################################################
# 1. DATA LOADING & COMBINING ACROSS TIMEPOINTS
##########################################################################################

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
cat("  Combined abundance table:", nrow(abund_combined), "features x", ncol(abund_combined) - 1, "sample-timepoints\n")

group_combined <- bind_rows(lapply(group_files, read_excel))
cat("  Combined group table:", nrow(group_combined), "sample-timepoints\n")

# Pre-clean constant (all-zero/all-NA) features and samples ourselves: LipidSigR's
# as_summarized_experiment() removes these internally but its rowname bookkeeping gets
# out of sync when BOTH axes have constant entries at once, causing a hard crash.
mat <- as.matrix(abund_combined[, -1])
rownames(mat) <- abund_combined$feature
feat_const <- apply(mat, 1, function(x) length(unique(x)) <= 1)
samp_const <- apply(mat, 2, function(x) length(unique(x)) <= 1)
if (any(feat_const)) cat("  Dropping", sum(feat_const), "constant-value feature(s):", paste(abund_combined$feature[feat_const], collapse = ", "), "\n")
if (any(samp_const)) cat("  Dropping", sum(samp_const), "constant-value sample-timepoint(s):", paste(colnames(mat)[samp_const], collapse = ", "), "\n")

abund_clean <- abund_combined %>% dplyr::filter(!feat_const)
keep_samples <- colnames(mat)[!samp_const]
abund_clean <- abund_clean[, c("feature", keep_samples)]
group_clean <- group_combined %>% dplyr::filter(sample_name %in% keep_samples)


##########################################################################################
# 2. LIPID ANNOTATION (rgoslin) + SUMMARIZEDEXPERIMENT + PROCESSING
##########################################################################################

cat("\n=== 2. LIPID ANNOTATION & PROCESSING ===\n")

parse_lipid <- rgoslin::parseLipidNames(lipidNames = abund_clean$feature)
recognized_lipids <- parse_lipid$Original.Name[which(parse_lipid$Grammar != "NOT_PARSEABLE")]
abundance_filtered <- abund_clean %>% dplyr::filter(feature %in% recognized_lipids)
goslin_annotation <- parse_lipid %>% dplyr::filter(Original.Name %in% recognized_lipids)
cat("  Recognized lipids:", length(recognized_lipids), "/", nrow(abund_clean), "\n")

# se_type = "profiling" is used deliberately: our own severity/timepoint grouping is
# handled outside LipidSigR (via the mixed model below), so we don't need the package's
# de_two/de_multiple group-comparison machinery here.
se <- as_summarized_experiment(
  abundance_filtered, goslin_annotation,
  group_info = group_clean, se_type = "profiling", paired_sample = NULL
)

processed_se <- data_process(
  se,
  exclude_missing = TRUE, exclude_missing_pct = 70,
  replace_na_method = "min", replace_na_method_ref = 0.5,
  normalization = "Percentage", transform = "log10"
)

data_process_plots <- plot_data_process(se, processed_se)
png(file.path(base_dir, "00.Data_quality/BoxPlot_before_process.png"), width = 1600, height = 1000, res = 150)
print(data_process_plots$static_boxPlot_before); dev.off()
png(file.path(base_dir, "00.Data_quality/BoxPlot_after_process.png"), width = 1600, height = 1000, res = 150)
print(data_process_plots$static_boxPlot_after); dev.off()

processed_abund <- processed_se@metadata$processed_abund
cat("  Processed (normalized, log10) lipids retained:", nrow(processed_abund), "\n")


##########################################################################################
# 3. BUILD PATIENT/TIMEPOINT/SEVERITY LONG-FORMAT TABLE
##########################################################################################

cat("\n=== 3. BUILDING LONG-FORMAT LONGITUDINAL TABLE ===\n")

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

cat("  Patients:", n_distinct(sample_meta$patient_id), "\n")
cat("  Timepoints per patient:\n")
print(table(table(sample_meta$patient_id)))

# Day of fever before hospitalization/sampling (per-patient, constant across timepoints).
# All mild/severe patients here were hospitalized -- this measures how far into their own
# illness course each patient already was at D0, which D0-as-study-day doesn't capture.
# Used below as a covariate so severity/time effects aren't confounded with this.
day_of_fever <- read_tsv("data/sick_patients_day_of_fever.tsv", show_col_types = FALSE)
sample_meta <- sample_meta %>% left_join(day_of_fever, by = c("patient_id" = "patient"))
n_missing_fever <- sum(is.na(sample_meta$day_of_fever))
if (n_missing_fever > 0) {
  cat("  Warning:", n_missing_fever, "sample-timepoint(s) missing day_of_fever (patient not in sick_patients_day_of_fever.tsv)\n")
}
fever_by_severity <- sample_meta %>% distinct(patient_id, severity, day_of_fever) %>% group_by(severity) %>%
  summarise(n = n(), mean_day_of_fever = mean(day_of_fever, na.rm = TRUE), sd = sd(day_of_fever, na.rm = TRUE))
cat("  Day of fever before hospitalization, by severity:\n")
print(fever_by_severity)

long_df <- processed_abund %>%
  pivot_longer(-feature, names_to = "sample_name", values_to = "value") %>%
  inner_join(sample_meta, by = "sample_name") %>%
  mutate(timepoint = factor(timepoint, levels = c("D0", "D3", "D10", "D60")))

write_tsv(long_df, file.path(base_dir, "02.Mixed_models/longitudinal_long_format.tsv"))


##########################################################################################
# 4. PER-LIPID LINEAR MIXED MODELS (timepoint * severity, patient random intercept)
##########################################################################################

cat("\n=== 4. FITTING PER-LIPID MIXED MODELS ===\n")

features <- unique(long_df$feature)
cat("  Fitting models for", length(features), "lipids...\n")

# Generic per-feature LMM fitter, reused below for class- and category-level aggregation
# (Section 7) so the exact same modeling logic applies at every level. day_of_fever is kept
# as a covariate in every one of the 4 models (never the term being tested), so it controls
# for baseline illness-course heterogeneity without changing what each LRT comparison isolates.
fit_one_feature <- function(f, df) {
  d <- df %>% filter(feature == f)
  out <- tryCatch({
    suppressMessages(withCallingHandlers({
      m_full   <- lmer(value ~ timepoint * severity + day_of_fever + (1 | patient_id), data = d, REML = FALSE)
      m_noint  <- lmer(value ~ timepoint + severity   + day_of_fever + (1 | patient_id), data = d, REML = FALSE)
      m_notime <- lmer(value ~ severity                + day_of_fever + (1 | patient_id), data = d, REML = FALSE)
      m_nosev  <- lmer(value ~ timepoint                + day_of_fever + (1 | patient_id), data = d, REML = FALSE)

      p_interaction <- anova(m_full, m_noint)[["Pr(>Chisq)"]][2]
      p_time        <- anova(m_noint, m_notime)[["Pr(>Chisq)"]][2]
      p_severity    <- anova(m_noint, m_nosev)[["Pr(>Chisq)"]][2]

      data.frame(
        feature = f, n_obs = nrow(d), n_patients = n_distinct(d$patient_id),
        p_time = p_time, p_severity = p_severity, p_interaction = p_interaction,
        singular_fit = isSingular(m_full), error = NA_character_
      )
    }, warning = function(w) invokeRestart("muffleWarning")))
  }, error = function(e) {
    data.frame(
      feature = f, n_obs = nrow(d), n_patients = n_distinct(d$patient_id),
      p_time = NA_real_, p_severity = NA_real_, p_interaction = NA_real_,
      singular_fit = NA, error = conditionMessage(e)
    )
  })
  out
}

lmm_results <- bind_rows(lapply(features, fit_one_feature, df = long_df))

n_failed <- sum(!is.na(lmm_results$error))
if (n_failed > 0) cat("  Warning:", n_failed, "lipid(s) failed to fit:", paste(lmm_results$feature[!is.na(lmm_results$error)], collapse = ", "), "\n")
cat("  Singular fits (near-zero patient variance, interpret with caution):", sum(lmm_results$singular_fit, na.rm = TRUE), "\n")

lmm_results <- lmm_results %>%
  mutate(
    fdr_time = p.adjust(p_time, method = "BH"),
    fdr_severity = p.adjust(p_severity, method = "BH"),
    fdr_interaction = p.adjust(p_interaction, method = "BH")
  ) %>%
  arrange(fdr_time)

write_tsv(lmm_results, file.path(base_dir, "02.Mixed_models/LMM_all_results.tsv"))

sig_time <- lmm_results %>% filter(fdr_time < 0.05) %>% arrange(fdr_time)
sig_interaction <- lmm_results %>% filter(fdr_interaction < 0.05) %>% arrange(fdr_interaction)
write_tsv(sig_time, file.path(base_dir, "02.Mixed_models/LMM_significant_time_effect.tsv"))
write_tsv(sig_interaction, file.path(base_dir, "02.Mixed_models/LMM_significant_interaction.tsv"))

cat("  Lipids with significant TIME effect (FDR<0.05):", nrow(sig_time), "\n")
cat("  Lipids with significant TIMEPOINT:SEVERITY interaction (FDR<0.05):", nrow(sig_interaction), "\n")


##########################################################################################
# 4b. POST-HOC PAIRWISE CONTRASTS (which timepoint drives the divergence?)
##########################################################################################
# The LRTs above only say a lipid's trajectory differs by severity (or changes over time),
# not WHICH timepoint(s) it differs at. For the lipids already flagged as significant, fit
# m_full again (cheap: only a handful of lipids, not all 119) and use emmeans to contrast
# each non-baseline timepoint against D0, separately within mild and within severe. This
# is the post-hoc step proposed but not implemented in the report's "Possible next steps".

cat("\n=== 4b. POST-HOC PAIRWISE CONTRASTS (emmeans, D0 as reference) ===\n")

fit_posthoc_contrasts <- function(f, df) {
  d <- df %>% filter(feature == f)
  tryCatch({
    suppressMessages(withCallingHandlers({
      m_full <- lmer(value ~ timepoint * severity + day_of_fever + (1 | patient_id), data = d, REML = FALSE)
      emm <- emmeans(m_full, ~ timepoint | severity)
      ct <- as.data.frame(contrast(emm, method = "trt.vs.ctrl", ref = 1))
      ct %>%
        transmute(
          feature = f, severity = severity, contrast = contrast,
          estimate = estimate, se = SE, p_value = p.value
        )
    }, warning = function(w) invokeRestart("muffleWarning")))
  }, error = function(e) {
    data.frame(feature = f, severity = NA_character_, contrast = NA_character_,
      estimate = NA_real_, se = NA_real_, p_value = NA_real_)
  })
}

posthoc_features <- union(sig_time$feature, sig_interaction$feature)
cat("  Fitting post-hoc contrasts for", length(posthoc_features), "already-significant lipid(s)...\n")

posthoc_results <- bind_rows(lapply(posthoc_features, fit_posthoc_contrasts, df = long_df)) %>%
  mutate(fdr_posthoc = p.adjust(p_value, method = "BH")) %>%
  arrange(fdr_posthoc)

write_tsv(posthoc_results, file.path(base_dir, "02.Mixed_models/LMM_posthoc_contrasts.tsv"))
cat("  Saved post-hoc contrasts (", nrow(posthoc_results), "rows: ", length(posthoc_features),
    "lipids x 2 severities x 3 timepoints)\n")


##########################################################################################
# 5. TRAJECTORY (SPAGHETTI) PLOTS FOR TOP HITS
##########################################################################################

cat("\n=== 5. TRAJECTORY PLOTS ===\n")

plot_trajectories <- function(feats, title, outfile, data = long_df) {
  if (length(feats) == 0) return(invisible(NULL))
  d <- data %>% filter(feature %in% feats) %>% mutate(feature = factor(feature, levels = feats))
  p <- ggplot(d, aes(x = timepoint, y = value, group = patient_id, color = severity)) +
    geom_line(alpha = 0.35) +
    geom_point(alpha = 0.5, size = 1) +
    stat_summary(aes(group = severity), fun = mean, geom = "line", linewidth = 1.3) +
    facet_wrap(~feature, scales = "free_y") +
    labs(title = title, x = "Timepoint", y = "log10(normalized abundance)", color = "Severity") +
    theme_bw()
  n_panels <- length(feats)
  ncol_panels <- min(4, n_panels)
  nrow_panels <- ceiling(n_panels / ncol_panels)
  png(outfile, width = 400 * ncol_panels, height = 350 * nrow_panels, res = 130)
  print(p)
  dev.off()
}

top_time <- head(sig_time$feature, 12)
top_interaction <- head(sig_interaction$feature, 12)

plot_trajectories(top_time, "Top lipids: significant TIME effect",
  file.path(base_dir, "02.Mixed_models/Trajectory_plots/Top_time_effect.png"))
plot_trajectories(top_interaction, "Top lipids: significant TIMEPOINT x SEVERITY interaction",
  file.path(base_dir, "02.Mixed_models/Trajectory_plots/Top_interaction_effect.png"))

cat("  Saved trajectory plots for top", length(top_time), "time-effect and", length(top_interaction), "interaction-effect lipids\n")

# The plots above are on the log10(% of total lipidome) scale used for modeling, which is
# hard to read directionally (values are negative whenever a lipid is <1% of the total).
# For interpretation, also show each patient's trajectory as log2 fold-change from their
# own D0 baseline: 0 = no change, positive = increase, negative = decrease. This does not
# change the statistics above, only how the same modeled values are displayed. Reused as-is
# for class/category level in Section 7.
make_log2fc <- function(df) {
  df %>%
    group_by(feature, patient_id) %>%
    mutate(baseline = value[timepoint == "D0"][1]) %>%
    ungroup() %>%
    mutate(log2fc = (value - baseline) / log10(2))
}
log2fc_df <- make_log2fc(long_df)

plot_fc_trajectories <- function(feats, title, outfile, data = log2fc_df) {
  if (length(feats) == 0) return(invisible(NULL))
  d <- data %>% filter(feature %in% feats) %>% mutate(feature = factor(feature, levels = feats))
  p <- ggplot(d, aes(x = timepoint, y = log2fc, group = patient_id, color = severity)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_line(alpha = 0.35) +
    geom_point(alpha = 0.5, size = 1) +
    stat_summary(aes(group = severity), fun = mean, geom = "line", linewidth = 1.3) +
    facet_wrap(~feature, scales = "free_y") +
    labs(title = title, x = "Timepoint", y = "log2 fold-change vs D0", color = "Severity") +
    theme_bw()
  n_panels <- length(feats)
  ncol_panels <- min(4, n_panels)
  nrow_panels <- ceiling(n_panels / ncol_panels)
  png(outfile, width = 400 * ncol_panels, height = 350 * nrow_panels, res = 130)
  print(p)
  dev.off()
}

plot_fc_trajectories(top_time, "Top lipids: significant TIME effect (fold-change vs D0)",
  file.path(base_dir, "02.Mixed_models/Trajectory_plots/Top_time_effect_FCfromD0.png"))
plot_fc_trajectories(top_interaction, "Top lipids: significant TIMEPOINT x SEVERITY interaction (fold-change vs D0)",
  file.path(base_dir, "02.Mixed_models/Trajectory_plots/Top_interaction_effect_FCfromD0.png"))

write_tsv(log2fc_df, file.path(base_dir, "02.Mixed_models/longitudinal_log2FC_from_D0.tsv"))

cat("  Saved fold-change-from-D0 trajectory plots\n")

# Effect-size summary: for each lipid, the largest |group-mean log2FC vs D0| seen across the
# 3 non-baseline timepoints x 2 severities. Group means (not per-patient values) are used so
# one outlier patient can't inflate the number. This gives a magnitude to sort/filter on,
# alongside the p-values already in LMM_all_results.tsv (significance alone isn't enough to
# prioritize biomarker candidates).
effect_size_summary <- log2fc_df %>%
  filter(timepoint != "D0") %>%
  group_by(feature, severity, timepoint) %>%
  summarise(mean_log2fc = mean(log2fc, na.rm = TRUE), .groups = "drop") %>%
  group_by(feature) %>%
  summarise(max_abs_group_log2fc = max(abs(mean_log2fc), na.rm = TRUE), .groups = "drop")

lmm_results <- lmm_results %>%
  left_join(effect_size_summary, by = "feature") %>%
  arrange(fdr_time)
write_tsv(lmm_results, file.path(base_dir, "02.Mixed_models/LMM_all_results.tsv"))
cat("  Added max_abs_group_log2fc effect-size column to LMM_all_results.tsv\n")


##########################################################################################
# 6. PCA TRAJECTORY PLOT (each patient's path through lipidomic space over time)
##########################################################################################

cat("\n=== 6. PCA TRAJECTORY ===\n")

result_pca <- dr_pca(
  processed_se, scaling = TRUE, centering = TRUE,
  clustering = "kmeans", cluster_num = 2,
  kmedoids_metric = NULL, distfun = NULL, hclustfun = NULL, eps = NULL, minPts = NULL,
  feature_contrib_pc = c(1, 2), plot_topN = 10
)

pca_scores <- result_pca$pca_rotated_data %>%
  rename(sample_name = 1) %>%
  inner_join(sample_meta, by = "sample_name") %>%
  mutate(timepoint = factor(timepoint, levels = c("D0", "D3", "D10", "D60")))

write_tsv(pca_scores, file.path(base_dir, "01.Trajectories/PCA/table_pca_scores_with_metadata.tsv"))

p_traj <- ggplot(pca_scores, aes(x = PC1, y = PC2, group = patient_id, color = severity)) +
  geom_path(arrow = arrow(length = unit(0.12, "cm")), alpha = 0.5) +
  geom_point(aes(shape = timepoint), size = 2) +
  facet_wrap(~severity) +
  labs(title = "Patient trajectories through PCA space (D0 -> D3 -> D10 -> D60)") +
  theme_bw()

png(file.path(base_dir, "01.Trajectories/PCA/Patient_trajectories_PCA.png"), width = 2000, height = 1000, res = 150)
print(p_traj)
dev.off()

cat("  Saved PCA trajectory plot\n")


##########################################################################################
# 7. LIPID CLASS / CATEGORY (FAMILY) LEVEL AGGREGATION
##########################################################################################
#
# Per-species trajectories (Section 4-5) are numerous and some are noisy/low-abundance.
# Aggregating to rgoslin's own lipid classification gives fewer, more interpretable
# trajectories at two levels:
#   - "class"    = rgoslin Lipid.Maps.Main.Class    (e.g. PC, SM, TG, DG, CE, LPC, FA, CAR)
#   - "category" = rgoslin Lipid.Maps.Category       (the broad families: Fatty Acyls,
#                  Glycerolipids, Glycerophospholipids, Sphingolipids, Sterol Lipids)
# Species-level values are on the log10(percentage-of-total) scale, so they cannot simply
# be summed (sum of logs != log of sum). Values are converted back to linear percentage
# (10^value), summed within each class/category per sample, then re-log10-transformed.
# The same paired LMM (Section 4) is refit on these aggregated series, since the same
# pseudoreplication argument against plain ANOVA applies here too.

cat("\n=== 7. CLASS / CATEGORY-LEVEL AGGREGATION ===\n")

category_labels <- c(
  FA = "Fatty Acyls", GL = "Glycerolipids", GP = "Glycerophospholipids",
  SP = "Sphingolipids", ST = "Sterol Lipids"
)
feature_char <- goslin_annotation %>%
  transmute(
    feature = Original.Name,
    lipid_class = Lipid.Maps.Main.Class,
    lipid_category = unname(category_labels[Lipid.Maps.Category])
  ) %>%
  filter(feature %in% processed_abund$feature) %>%
  distinct(feature, .keep_all = TRUE)

linear_abund <- processed_abund
linear_abund[-1] <- 10^linear_abund[-1]

aggregate_char <- function(char_col) {
  linear_abund %>%
    inner_join(feature_char %>% select(feature, char = all_of(char_col)), by = "feature") %>%
    select(-feature) %>%
    group_by(char) %>%
    summarise(across(everything(), sum), .groups = "drop") %>%
    rename(feature = char) %>%
    mutate(across(-feature, ~ log10(.x)))
}

class_abund <- aggregate_char("lipid_class")
category_abund <- aggregate_char("lipid_category")
cat("  Aggregated", nrow(processed_abund), "species into", nrow(class_abund), "classes and", nrow(category_abund), "categories (families)\n")

class_long_df <- class_abund %>%
  pivot_longer(-feature, names_to = "sample_name", values_to = "value") %>%
  inner_join(sample_meta, by = "sample_name") %>%
  mutate(timepoint = factor(timepoint, levels = c("D0", "D3", "D10", "D60")))
category_long_df <- category_abund %>%
  pivot_longer(-feature, names_to = "sample_name", values_to = "value") %>%
  inner_join(sample_meta, by = "sample_name") %>%
  mutate(timepoint = factor(timepoint, levels = c("D0", "D3", "D10", "D60")))

write_tsv(class_long_df, file.path(base_dir, "03.Aggregated/Class/class_long_format.tsv"))
write_tsv(category_long_df, file.path(base_dir, "03.Aggregated/Category/category_long_format.tsv"))

run_aggregated_lmm <- function(df, out_dir, label) {
  feats <- unique(df$feature)
  res <- bind_rows(lapply(feats, fit_one_feature, df = df)) %>%
    mutate(
      fdr_time = p.adjust(p_time, method = "BH"),
      fdr_severity = p.adjust(p_severity, method = "BH"),
      fdr_interaction = p.adjust(p_interaction, method = "BH")
    ) %>%
    arrange(fdr_time)
  write_tsv(res, file.path(out_dir, paste0("LMM_", label, "_results.tsv")))
  cat("  [", label, "] significant TIME effect (FDR<0.05):", sum(res$fdr_time < 0.05, na.rm = TRUE),
      "/", nrow(res), " | significant interaction:", sum(res$fdr_interaction < 0.05, na.rm = TRUE), "\n")
  res
}

class_results <- run_aggregated_lmm(class_long_df, file.path(base_dir, "03.Aggregated/Class"), "class")
category_results <- run_aggregated_lmm(category_long_df, file.path(base_dir, "03.Aggregated/Category"), "category")

# Few enough classes/categories to plot all of them at once (not just top hits)
class_log2fc <- make_log2fc(class_long_df)
category_log2fc <- make_log2fc(category_long_df)

plot_trajectories(class_results$feature, "All lipid classes",
  file.path(base_dir, "03.Aggregated/Class/Trajectory_plots/All_classes.png"), data = class_long_df)
plot_fc_trajectories(class_results$feature, "All lipid classes (fold-change vs D0)",
  file.path(base_dir, "03.Aggregated/Class/Trajectory_plots/All_classes_FCfromD0.png"), data = class_log2fc)

plot_trajectories(category_results$feature, "All lipid categories (rgoslin families)",
  file.path(base_dir, "03.Aggregated/Category/Trajectory_plots/All_categories.png"), data = category_long_df)
plot_fc_trajectories(category_results$feature, "All lipid categories (fold-change vs D0)",
  file.path(base_dir, "03.Aggregated/Category/Trajectory_plots/All_categories_FCfromD0.png"), data = category_log2fc)

cat("  Saved class- and category-level trajectory plots\n")


##########################################################################################
# 8. LIGHTWEIGHT CLASS/CATEGORY OVER-REPRESENTATION CHECK
##########################################################################################
# LipidSigR's own enrichment_lsea()/enrichment_ora() require a deSp_se object built by
# deSp_multiGroup()/deSp_twoGroup(), which this pipeline deliberately doesn't use (see
# header comment - avoiding pseudoreplication). Rather than force our LMM results into that
# API, run a simple Fisher's exact test per class/category: is a given significant-feature
# set (e.g. the interaction-significant "severity-divergence" lipids) over-represented in
# any lipid family, relative to all tested lipids as background? A family-level signal is
# more biologically credible than scattered single-lipid hits.

cat("\n=== 8. CLASS/CATEGORY OVER-REPRESENTATION CHECK ===\n")

run_ora <- function(sig_features, background_features, char_df, char_col, label) {
  bg <- char_df %>% filter(feature %in% background_features)
  n_total <- nrow(bg)
  bg %>%
    distinct(.data[[char_col]]) %>%
    pull(1) %>%
    lapply(function(grp) {
      in_group <- bg[[char_col]] == grp
      in_set <- bg$feature %in% sig_features
      tbl <- table(in_set, in_group)
      ft <- fisher.test(tbl)
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

all_species <- unique(long_df$feature)
ora_results <- bind_rows(
  run_ora(sig_interaction$feature, all_species, feature_char, "lipid_class", "class (interaction hits)"),
  run_ora(sig_interaction$feature, all_species, feature_char, "lipid_category", "category (interaction hits)"),
  run_ora(sig_time$feature, all_species, feature_char, "lipid_class", "class (time hits)"),
  run_ora(sig_time$feature, all_species, feature_char, "lipid_category", "category (time hits)")
)

write_tsv(ora_results, file.path(base_dir, "03.Aggregated/Enrichment_interaction_hits.tsv"))
cat("  Saved class/category over-representation results (", nrow(ora_results), "group tests)\n")


##########################################################################################
# DONE
##########################################################################################

cat("\n========================================\n")
cat("LONGITUDINAL ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("Output directory:", base_dir, "\n")
cat("  Species-level results:  02.Mixed_models/\n")
cat("  Post-hoc contrasts:     02.Mixed_models/LMM_posthoc_contrasts.tsv\n")
cat("  Class-level results:    03.Aggregated/Class/\n")
cat("  Category-level results: 03.Aggregated/Category/\n")
cat("  Over-representation:   03.Aggregated/Enrichment_interaction_hits.tsv\n")
