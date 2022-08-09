# Short script to generate a one document report from all pdfs in the docs 
# folder
# Jakob J. Assmann j.assmann@bio.au.dk 2 March 2022

# Dependencies
library(rmarkdown)
library(webshot)
library(pdftools)

# Set base names and order of files to include
report_files <- paste0("docs/",
                        c("index.pdf",
                          "workflow.pdf",
                          "training_annotations.pdf",
                          "data_overview.pdf",
                          "focal_var_selection.pdf",
                          "gbm_models_performance.pdf",
                          "ranger_models_performance.pdf",
                          "summary_stats.pdf"
                        ))

# Generate pdf of index.md
render("docs/index.md", html_document())

# Check whether other files whether an RMD exist, rerender and convert to PDF
lapply(report_files, function(x){
  x <- gsub("(.*)\\..*", "\\1", x)
  cat("Processing:", x, "\n")
  # if(file.exists(paste0(x, ".Rmd"))) {
  #   cat("\tRMD file found re-generating...\n")
  #   #render(paste0(x, ".Rmd"), quiet = T)
  # }
  cat("\tChecking and removing old PDF..\n")
  if(file.exists(paste0(x, ".pdf"))) file.remove(paste0(x, ".pdf"))
  cat("\tGenerating new PDF...")
  webshot(paste0(x, ".html"),
          paste0(x, ".pdf"))
  cat("done.\n")
  return("OK")
})

# Combine to one report
pdf_combine(report_files,
             output = "docs/Assmann_et_al-DK_Forest_Quality_Report_v0.9.1.pdf")

# Remove intermediate pdf files
lapply(report_files, file.remove)
file.remove("docs/index.html")
