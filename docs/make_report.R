# Short script to generate a one document report from all pdfs in the docs 
# folder
# Jakob J. Assmann j.assmann@bio.au.dk 2 March 2022

# Dependencies
library(rmarkdown)
library(pdftools)

# Generate pdf of index.md
render("index.md", pdf_document())

# Set list of pdfs to combine
report_files <- c(
  "index.pdf",
  "data_overview.pdf",
  "focal_var_selection.pdf",
  "gbm_models_performance.pdf",
  "ranger_models_performance.pdf",
  "summary_stats.pdf"
)

# Combine PDFs using pdf tools
pdf_combined(report_files,
             output = "Assmann_et_al-DK_Forest_Quality_v0.1.0")