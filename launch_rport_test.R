#!/usr/bin/env /opt/R/4.3.1/bin/Rscript
# Save the entire workspace
suppressPackageStartupMessages(library(hipathia))
source("src/report.R")

output_folder <- "tmp/metabopathia_report5/"
load(file = file.path(output_folder,"workspace.RData"))
attach(opt)

servr::daemon_stop() # kill & close
# DAreport() to easily create a report
# Save and serve all results to browser
MTreport <- DAreport(met_results, metabo_pathways, path = output_folder, output_folder = "metabopathia_report", verbose = verbose)
status("100", paste0("HTML report created successfully in the ",file.path(codebase,MTreport)), output_folder)

serve_report(file.path(codebase,MTreport), port = servr::random_port(), browser = T, daemon = T)
message("Press Ctrl + C to stop serving the report...\n")
if(hipathia){
  HIreport <- DAreport(hi_results, pathways, path = output_folder, output_folder = "hipathia4comp_report",verbose = verbose)
  message("Press Ctrl + C to stop serving the report...\n")
  serve_report(file.path(codebase,HIreport),port = servr::random_port(), browser = T, daemon = F)
}

# Pause and wait for user input
if(!is_rstudio){
  cat("Press Enter to finish...")
  invisible(readLines("stdin", n=1))
}


