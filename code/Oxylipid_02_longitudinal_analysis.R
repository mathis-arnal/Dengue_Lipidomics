##########################################################################################
# OXYLIPID LONGITUDINAL ANALYSIS: D0 -> D3 -> D10 -> D60 (MILD vs SEVERE, PAIRED)
##########################################################################################
#
# GOAL: track how each patient's oxylipin profile evolves across all 4 timepoints and
# resolves toward its own D60 (convalescence) level, rather than comparing groups
# cross-sectionally at a single timepoint. This is the oxylipid analogue of
# 05_longitudinal_analysis.R, built on the same design: a per-lipid linear mixed model
# with patient as a random intercept, rather than treating each patient's 4 repeated
# measurements as independent observations.
#
# DATA SOURCE: data/lipidsig_datasets/Oxylipines_raw/Oxylipid_lipid_abundance_data_*.tsv +
# group_information_table_oxylipid_*.tsv, built by code/Oxylipid_lipidomics_pretreatment.R
# from the raw (pre-ratio) abundances in data/raw_data/oxylipines.xlsx. There is no
# healthy arm (healthy donors are single-draw, no D60 sample).
#
# WHY NOT LipidSigR's as_summarized_experiment() / data_process() / dr_pca()?
# Oxylipins aren't parseable by rgoslin (see Oxylipid_class_annotation.R header for the
# empirical crash this was tried and hit), so this script reimplements the same
# normalization (percentage + log10) and PCA directly in plain R and uses a manually
# curated class-annotation table (parent-PUFA family + compound class) in place of
# rgoslin's Lipid.Maps.Category / Lipid.Maps.Main.Class for the aggregation step.
#
# STATISTICAL MODEL (per lipid, identical in form to 05_longitudinal_analysis.R):
#   value ~ timepoint * severity + day_of_fever + (1 | patient_id)
# timepoint is a 4-level FACTOR (D0/D3/D10/D60 spacing is highly uneven: 0/3/10/60 days).
# day_of_fever (days of fever before hospitalization/sampling, from
# data/sick_patients_day_of_fever.tsv) is included as a covariate in every model so
# severity/time effects aren't confounded with how far into their illness course each
# patient already was at enrollment (D0). Models fit by ML (REML=FALSE) for LRT-based
# nested-model comparison:
#   - p_time        : does abundance change over time at all (averaged across groups)?
#   - p_severity    : do groups differ on average (across time)?
#   - p_interaction : does the trajectory/resolution shape differ between mild and severe?
# All p-values BH-FDR-adjusted across the 24 tested lipids.
##########################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(lme4)
  library(emmeans)
  library(ggplot2)
})

base_dir <- file.path(getwd(), "analysis", "Oxylipids", "Longitudinal")
data_dir <- "data/lipidsig_datasets/Oxylipines_raw"

dir.create(file.path(base_dir, "00.Data_quality"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "01.Trajectories/PCA"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "02.Mixed_models/Trajectory_plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "03.Aggregated/Family/Trajectory_plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "03.Aggregated/Class/Trajectory_plots"), recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("OXYLIPID LONGITUDINAL ANALYSIS (D0/D3/D10/D60)\n")
cat("Groups: Mild, Severe (paired, no healthy longitudinal data)\n")
cat("========================================\n")


##########################################################################################
# 1. DATA LOADING & COMBINING ACROSS TIMEPOINTS
##########################################################################################

cat("\n=== 1. DATA LOADING ===\n")

tp_suffix <- c(D0 = "D0", D3 = "D03", D10 = "D10", D60 = "D60")

abund_list <- lapply(tp_suffix, function(s) read_tsv(file.path(data_dir, sprintf("Oxylipid_lipid_abundance_data_%s.tsv", s)), show_col_types = FALSE))
stopifnot(all(vapply(abund_list, function(x) identical(x$feature, abund_list$D0$feature), logical(1))))

abund_combined <- abund_list$D0
for (nm in c("D3", "D10", "D60")) {
  abund_combined <- dplyr::full_join(abund_combined, abund_list[[nm]], by = "feature")
}
cat("  Combined abundance table:", nrow(abund_combined), "features x", ncol(abund_combined) - 1, "sample-timepoints\n")

group_list <- lapply(tp_suffix, function(s) read_tsv(file.path(data_dir, sprintf("group_information_table_oxylipid_%s.tsv", s)), show_col_types = FALSE))
group_combined <- bind_rows(group_list)
cat("  Combined group table:", nrow(group_combined), "sample-timepoints\n")

class_annotation <- read_tsv(file.path(data_dir, "oxylipin_class_annotation.tsv"), show_col_types = FALSE)

# Diagnostic only (not a required workaround here -- the as_summarized_experiment
# constant-row/column crash this mirrors in 05_longitudinal_analysis.R doesn't apply since
# we never build a SummarizedExperiment for oxylipids).
mat <- as.matrix(abund_combined[, -1])
feat_const <- apply(mat, 1, function(x) length(unique(x[!is.na(x)])) <= 1)
if (any(feat_const)) cat("  Note:", sum(feat_const), "constant-value feature(s):", paste(abund_combined$feature[feat_const], collapse = ", "), "\n")


##########################################################################################
# 2. NORMALIZATION (percentage + log10)
##########################################################################################

cat("\n=== 2. NORMALIZATION ===\n")

abund_mat <- as.matrix(abund_combined[, -1])
rownames(abund_mat) <- abund_combined$feature

# Drop features detected (non-zero) in fewer than exclude_missing_pct% of sample-timepoints
# -- matches data_process(exclude_missing=TRUE, exclude_missing_pct=70) elsewhere in this
# project. Needed for the expanded 47-compound panel: several resolvins/protectins/minor
# prostaglandins are undetectable in most samples, which would otherwise leave a constant
# (all-floor-value) row after normalize_log10's zero-replacement and break prcomp() later.
filter_low_prevalence <- function(m, exclude_missing_pct = 70) {
  detected_pct <- 100 * rowMeans(m > 0 & !is.na(m))
  keep <- !is.na(detected_pct) & detected_pct >= exclude_missing_pct
  if (any(!keep)) {
    cat("  Dropping", sum(!keep), "feature(s) detected in <", exclude_missing_pct, "% of sample-timepoints:",
        paste(rownames(m)[!keep], collapse = ", "), "\n")
  }
  m[keep, , drop = FALSE]
}
abund_mat <- filter_low_prevalence(abund_mat)

normalize_log10 <- function(m) {
  pct <- sweep(m, 2, colSums(m, na.rm = TRUE), FUN = "/") * 100
  pct_filled <- t(apply(pct, 1, function(row) {
    nonzero_min <- suppressWarnings(min(row[row > 0], na.rm = TRUE))
    if (!is.finite(nonzero_min)) nonzero_min <- 1e-6
    row[is.na(row) | row <= 0] <- nonzero_min * 0.5
    row
  }))
  dimnames(pct_filled) <- dimnames(pct)
  log10(pct_filled)
}

processed <- normalize_log10(abund_mat)
processed_abund <- as.data.frame(processed) %>% tibble::rownames_to_column("feature")

qc_before <- as.data.frame(abund_mat) %>% tibble::rownames_to_column("feature") %>%
  pivot_longer(-feature, names_to = "sample_name", values_to = "value") %>% mutate(stage = "before")
qc_after <- processed_abund %>%
  pivot_longer(-feature, names_to = "sample_name", values_to = "value") %>% mutate(stage = "after (log10 %)")

p_before <- ggplot(qc_before, aes(x = sample_name, y = value)) + geom_boxplot(outlier.size = 0.3) +
  labs(title = "Before normalization", x = NULL, y = "raw abundance") + theme_bw() +
  theme(axis.text.x = element_blank())
p_after <- ggplot(qc_after, aes(x = sample_name, y = value)) + geom_boxplot(outlier.size = 0.3) +
  labs(title = "After normalization", x = NULL, y = "log10(%)") + theme_bw() +
  theme(axis.text.x = element_blank())
png(file.path(base_dir, "00.Data_quality/BoxPlot_before_process.png"), width = 1800, height = 900, res = 150)
print(p_before); dev.off()
png(file.path(base_dir, "00.Data_quality/BoxPlot_after_process.png"), width = 1800, height = 900, res = 150)
print(p_after); dev.off()

cat("  Processed (normalized, log10):", nrow(processed_abund), "lipids x", ncol(processed_abund) - 1, "sample-timepoints\n")


##########################################################################################
# 3. BUILD PATIENT/TIMEPOINT/SEVERITY LONG-FORMAT TABLE
##########################################################################################

cat("\n=== 3. BUILDING LONG-FORMAT LONGITUDINAL TABLE ===\n")

sample_meta <- group_combined %>%
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

day_of_fever <- read_tsv("data/sick_patients_day_of_fever.tsv", show_col_types = FALSE)
sample_meta <- sample_meta %>% left_join(day_of_fever, by = c("patient_id" = "patient"))
n_missing_fever <- sum(is.na(sample_meta$day_of_fever))
if (n_missing_fever > 0) {
  cat("  Warning:", n_missing_fever, "sample-timepoint(s) missing day_of_fever (patient not in sick_patients_day_of_fever.tsv)\n")
}

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

if (length(posthoc_features) > 0) {
  posthoc_results <- bind_rows(lapply(posthoc_features, fit_posthoc_contrasts, df = long_df)) %>%
    mutate(fdr_posthoc = p.adjust(p_value, method = "BH")) %>%
    arrange(fdr_posthoc)
  write_tsv(posthoc_results, file.path(base_dir, "02.Mixed_models/LMM_posthoc_contrasts.tsv"))
  cat("  Saved post-hoc contrasts (", nrow(posthoc_results), "rows)\n")
} else {
  cat("  No lipids reached significance -- skipping post-hoc contrasts\n")
}


##########################################################################################
# 5. TRAJECTORY (SPAGHETTI) PLOTS
##########################################################################################

cat("\n=== 5. TRAJECTORY PLOTS ===\n")

plot_trajectories <- function(feats, title, outfile, data = long_df, y_lab = "log10(normalized abundance)", y_col = "value") {
  if (length(feats) == 0) return(invisible(NULL))
  d <- data %>% filter(feature %in% feats) %>% mutate(feature = factor(feature, levels = feats))
  p <- ggplot(d, aes(x = timepoint, y = .data[[y_col]], group = patient_id, color = severity)) +
    geom_line(alpha = 0.35) +
    geom_point(alpha = 0.5, size = 1) +
    stat_summary(aes(group = severity), fun = mean, geom = "line", linewidth = 1.3) +
    facet_wrap(~feature, scales = "free_y") +
    labs(title = title, x = "Timepoint", y = y_lab, color = "Severity") +
    theme_bw()
  n_panels <- length(feats)
  ncol_panels <- min(4, n_panels)
  nrow_panels <- ceiling(n_panels / ncol_panels)
  png(outfile, width = 400 * ncol_panels, height = 350 * nrow_panels, res = 130)
  print(p)
  dev.off()
}

# With only 24 lipids total, plot all of them (not just top hits) alongside the
# significant subset.
all_features_sorted <- lmm_results$feature
plot_trajectories(all_features_sorted, "All oxylipins",
  file.path(base_dir, "02.Mixed_models/Trajectory_plots/All_lipids.png"))
plot_trajectories(sig_time$feature, "Oxylipins: significant TIME effect",
  file.path(base_dir, "02.Mixed_models/Trajectory_plots/Significant_time_effect.png"))
plot_trajectories(sig_interaction$feature, "Oxylipins: significant TIMEPOINT x SEVERITY interaction",
  file.path(base_dir, "02.Mixed_models/Trajectory_plots/Significant_interaction_effect.png"))

cat("  Saved trajectory plots (all lipids, ", nrow(sig_time), "time-effect, ", nrow(sig_interaction), "interaction-effect)\n")

# log2 fold-change from each patient's own D0 baseline -- easier to read directionally
# than the log10(% of total) scale used for modeling. Reused for family/class level below.
make_log2fc <- function(df) {
  df %>%
    group_by(feature, patient_id) %>%
    mutate(baseline = value[timepoint == "D0"][1]) %>%
    ungroup() %>%
    mutate(log2fc = (value - baseline) / log10(2))
}
log2fc_df <- make_log2fc(long_df)

plot_trajectories(all_features_sorted, "All oxylipins (fold-change vs D0)",
  file.path(base_dir, "02.Mixed_models/Trajectory_plots/All_lipids_FCfromD0.png"),
  data = log2fc_df, y_lab = "log2 fold-change vs D0", y_col = "log2fc")

write_tsv(log2fc_df, file.path(base_dir, "02.Mixed_models/longitudinal_log2FC_from_D0.tsv"))
cat("  Saved fold-change-from-D0 trajectory plot\n")

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


##########################################################################################
# 6. PCA TRAJECTORY PLOT
##########################################################################################

cat("\n=== 6. PCA TRAJECTORY ===\n")

pca <- prcomp(t(processed), scale. = TRUE, center = TRUE)
var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_scores <- as.data.frame(pca$x[, 1:2]) %>%
  tibble::rownames_to_column("sample_name") %>%
  inner_join(sample_meta, by = "sample_name") %>%
  mutate(timepoint = factor(timepoint, levels = c("D0", "D3", "D10", "D60")))

write_tsv(pca_scores, file.path(base_dir, "01.Trajectories/PCA/table_pca_scores_with_metadata.tsv"))

p_traj <- ggplot(pca_scores, aes(x = PC1, y = PC2, group = patient_id, color = severity)) +
  geom_path(arrow = arrow(length = unit(0.12, "cm")), alpha = 0.5) +
  geom_point(aes(shape = timepoint), size = 2) +
  facet_wrap(~severity) +
  labs(title = "Patient trajectories through PCA space (D0 -> D3 -> D10 -> D60)",
       x = sprintf("PC1 (%.1f%%)", var_explained[1]), y = sprintf("PC2 (%.1f%%)", var_explained[2])) +
  theme_bw()

png(file.path(base_dir, "01.Trajectories/PCA/Patient_trajectories_PCA.png"), width = 2000, height = 1000, res = 150)
print(p_traj)
dev.off()

cat("  Saved PCA trajectory plot\n")


##########################################################################################
# 7. PARENT-PUFA-FAMILY / CLASS LEVEL AGGREGATION
##########################################################################################
# Analogous to 05_longitudinal_analysis.R's rgoslin class/category aggregation, but using
# the manual oxylipin annotation (Oxylipid_class_annotation.R) since rgoslin doesn't apply:
#   - "class"  = compound class (HETE, HODE, PG, TX, ...)
#   - "family" = parent PUFA the oxylipin derives from (LA/AA/ALA/EPA/DHA-derived)
# Species-level values are log10(percentage), so they're converted back to linear percentage,
# summed within each class/family per sample, then re-log10-transformed before refitting the
# same LMM.

cat("\n=== 7. FAMILY / CLASS-LEVEL AGGREGATION ===\n")

linear_abund <- processed_abund
linear_abund[-1] <- 10^linear_abund[-1]

aggregate_char <- function(char_col) {
  linear_abund %>%
    inner_join(class_annotation %>% select(feature, char = all_of(char_col)), by = "feature") %>%
    select(-feature) %>%
    group_by(char) %>%
    summarise(across(everything(), sum), .groups = "drop") %>%
    rename(feature = char) %>%
    mutate(across(-feature, ~ log10(.x)))
}

family_abund <- aggregate_char("family")
class_abund <- aggregate_char("class")
cat("  Aggregated", nrow(processed_abund), "species into", nrow(family_abund), "families and", nrow(class_abund), "classes\n")

family_long_df <- family_abund %>%
  pivot_longer(-feature, names_to = "sample_name", values_to = "value") %>%
  inner_join(sample_meta, by = "sample_name") %>%
  mutate(timepoint = factor(timepoint, levels = c("D0", "D3", "D10", "D60")))
class_long_df <- class_abund %>%
  pivot_longer(-feature, names_to = "sample_name", values_to = "value") %>%
  inner_join(sample_meta, by = "sample_name") %>%
  mutate(timepoint = factor(timepoint, levels = c("D0", "D3", "D10", "D60")))

write_tsv(family_long_df, file.path(base_dir, "03.Aggregated/Family/family_long_format.tsv"))
write_tsv(class_long_df, file.path(base_dir, "03.Aggregated/Class/class_long_format.tsv"))

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

family_results <- run_aggregated_lmm(family_long_df, file.path(base_dir, "03.Aggregated/Family"), "family")
class_results <- run_aggregated_lmm(class_long_df, file.path(base_dir, "03.Aggregated/Class"), "class")

family_log2fc <- make_log2fc(family_long_df)
class_log2fc <- make_log2fc(class_long_df)

plot_trajectories(family_results$feature, "All parent-PUFA families",
  file.path(base_dir, "03.Aggregated/Family/Trajectory_plots/All_families.png"), data = family_long_df)
plot_trajectories(family_results$feature, "All parent-PUFA families (fold-change vs D0)",
  file.path(base_dir, "03.Aggregated/Family/Trajectory_plots/All_families_FCfromD0.png"),
  data = family_log2fc, y_lab = "log2 fold-change vs D0", y_col = "log2fc")

plot_trajectories(class_results$feature, "All oxylipin classes",
  file.path(base_dir, "03.Aggregated/Class/Trajectory_plots/All_classes.png"), data = class_long_df)
plot_trajectories(class_results$feature, "All oxylipin classes (fold-change vs D0)",
  file.path(base_dir, "03.Aggregated/Class/Trajectory_plots/All_classes_FCfromD0.png"),
  data = class_log2fc, y_lab = "log2 fold-change vs D0", y_col = "log2fc")

cat("  Saved family- and class-level trajectory plots\n")


##########################################################################################
# 8. LIGHTWEIGHT CLASS/FAMILY OVER-REPRESENTATION CHECK (Fisher's exact test)
##########################################################################################
# Same pattern as 05_longitudinal_analysis.R Section 8: is a significant-lipid set
# over-represented in any class/family, relative to all 24 tested lipids as background?

cat("\n=== 8. CLASS/FAMILY OVER-REPRESENTATION CHECK ===\n")

run_ora <- function(sig_features, background_features, char_df, char_col, label) {
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

all_species <- unique(long_df$feature)
ora_results <- bind_rows(
  run_ora(sig_interaction$feature, all_species, class_annotation, "class", "class (interaction hits)"),
  run_ora(sig_interaction$feature, all_species, class_annotation, "family", "family (interaction hits)"),
  run_ora(sig_time$feature, all_species, class_annotation, "class", "class (time hits)"),
  run_ora(sig_time$feature, all_species, class_annotation, "family", "family (time hits)")
)

write_tsv(ora_results, file.path(base_dir, "03.Aggregated/Enrichment_significant_hits.tsv"))
cat("  Saved class/family over-representation results (", nrow(ora_results), "group tests)\n")


##########################################################################################
# DONE
##########################################################################################

cat("\n========================================\n")
cat("OXYLIPID LONGITUDINAL ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("Output directory:", base_dir, "\n")
cat("  Species-level results:  02.Mixed_models/\n")
cat("  Post-hoc contrasts:     02.Mixed_models/LMM_posthoc_contrasts.tsv\n")
cat("  Family-level results:   03.Aggregated/Family/\n")
cat("  Class-level results:    03.Aggregated/Class/\n")
cat("  Over-representation:    03.Aggregated/Enrichment_significant_hits.tsv\n")
