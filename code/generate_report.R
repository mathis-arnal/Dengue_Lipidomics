##########################################################################################
# GENERATE HTML REPORT FROM RMD
##########################################################################################
#
# This script renders the R Markdown report to HTML format
#
##########################################################################################

# Load required library
library(rmarkdown)

# Set working directory to analysis folder (where images are located)
setwd("/Users/loictalignani/research/project/lipidomics/analysis")

# Render the report
rmarkdown::render(
  input = "Lipidomics_Analysis_Report.Rmd",
  output_format = "html_document",
  output_file = "Lipidomics_Analysis_Report.html",
  output_dir = "/Users/loictalignani/research/project/lipidomics/analysis"
)

cat("\n========================================\n")
cat("Report generated successfully!\n")
cat("Location: /Users/loictalignani/research/project/lipidomics/analysis/Lipidomics_Analysis_Report.html\n")
cat("========================================\n")
