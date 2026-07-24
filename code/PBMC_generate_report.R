##########################################################################################
# GENERATE HTML REPORT FROM RMD
##########################################################################################
#
# This script renders the R Markdown report to HTML format
#
##########################################################################################


source

# Load required library
library(rmarkdown)

setwd("E:/Dengue_lipidomics")

source("code/PBMC_01_analysis_three_groups_all.R")
source("code/PBMC_01_analysis_three_groups.R")
source("code/PBMC_02_posthoc_healthy_vs_mild.R")
source("code/PBMC_03_posthoc_healthy_vs_severe.R")
source("code/PBMC_04_posthoc_mild_vs_severe.R")
source("code/PBMC_05_total_abundance_check.R")

# Render the report
rmarkdown::render(
  input = "analysis/PBMC_Lipidomics_Analysis_Report.Rmd",
  output_format = "html_document",
  output_file = "PBMC_Lipidomics_Analysis_Report.html",
  output_dir = "analysis"
)

cat("\n========================================\n")
cat("Report generated successfully!\n")
cat("Location: E:/Dengue_lipidomics/analysis/PBMC_Lipidomics_Analysis_Report.html\n")
cat("========================================\n")
