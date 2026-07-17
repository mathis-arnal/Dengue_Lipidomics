##########################################################################################
# PBMC LIPIDOMICS: TOTAL ABUNDANCE BY SEVERITY (COMPLEMENTARY CHECK)
##########################################################################################
#
# LipidSigR's data_process(normalization = "Percentage") in PBMC_01-04 rescales every
# sample to its own total signal, which removes any between-sample difference in total
# lipid abundance before the ANOVA/t-tests run. That total-abundance axis is where a
# dengue-associated hypolipidemia signal would actually live, so it's checked here
# directly on the pre-normalization (raw) PBMC abundance table, independent of the
# species-level LipidSigR pipeline.
#
# Run from the repository root (relative paths, no setwd()).
##########################################################################################

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

base_dir <- "analysis/PBMC/Three_groups/D0/04.Total_abundance_check"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

group_info <- read_tsv("data/PBMC/group_information_table_healthy_vs_sick_patients_D0_all.tsv")
abundance <- read_tsv("data/PBMC/PBMC_Lipid_abundance_data_D0.tsv")

matched_samples <- intersect(group_info$sample_name, colnames(abundance))
abundance <- abundance %>% dplyr::select(feature, all_of(matched_samples))
group_info <- group_info %>% dplyr::filter(sample_name %in% matched_samples)

cat("Matched samples:", length(matched_samples), "\n")
cat("  Healthy:", sum(group_info$group == "healthy"), "\n")
cat("  Mild:", sum(group_info$group == "mild"), "\n")
cat("  Severe:", sum(group_info$group == "severe"), "\n")

totals <- abundance %>%
  dplyr::select(-feature) %>%
  summarise(across(everything(), \(x) sum(x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "sample_name", values_to = "total_abundance") %>%
  left_join(group_info %>% dplyr::select(sample_name, group), by = "sample_name") %>%
  mutate(group = factor(group, levels = c("healthy", "mild", "severe")))

write_tsv(totals, file.path(base_dir, "total_abundance_by_sample.tsv"))

summary_tbl <- totals %>%
  group_by(group) %>%
  summarise(
    n = n(),
    median = median(total_abundance),
    mean = mean(total_abundance),
    sd = sd(total_abundance),
    .groups = "drop"
  )
write_tsv(summary_tbl, file.path(base_dir, "total_abundance_summary_by_group.tsv"))

kw <- kruskal.test(total_abundance ~ group, data = totals)
kw_df <- data.frame(
  test = "Kruskal-Wallis (3 groups)",
  statistic = unname(kw$statistic),
  df = unname(kw$parameter),
  p.value = kw$p.value
)
write_tsv(kw_df, file.path(base_dir, "kruskal_wallis_result.tsv"))

pairwise <- pairwise.wilcox.test(totals$total_abundance, totals$group, p.adjust.method = "BH")
pw_df <- as.data.frame(as.table(pairwise$p.value)) %>%
  dplyr::filter(!is.na(Freq)) %>%
  dplyr::rename(group1 = Var1, group2 = Var2, p.adj = Freq)
write_tsv(pw_df, file.path(base_dir, "pairwise_wilcoxon_results.tsv"))

cat("\nKruskal-Wallis (healthy vs mild vs severe): p =", signif(kw$p.value, 3), "\n")
cat("\nPairwise Wilcoxon (BH-adjusted):\n")
print(pw_df)

p <- ggplot(totals, aes(x = group, y = total_abundance, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  scale_y_log10() +
  labs(
    title = "Total raw lipid abundance by severity group (PBMC, D0)",
    subtitle = sprintf("Kruskal-Wallis p = %.3f", kw$p.value),
    x = NULL,
    y = "Total abundance per sample (sum of raw feature values, log10 scale)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave(file.path(base_dir, "Total_abundance_by_severity_boxplot.png"), p, width = 6, height = 5, dpi = 150)

cat("\nOutputs written to:", base_dir, "\n")
