#!/usr/bin/env Rscript
# Install required packages for haplotype clustering tool.
# Run once: Rscript install.R

pkgs <- c(
  "vcfR",
  "ggplot2",
  "dplyr",
  "tidyr",
  "tibble",
  "pheatmap",
  "RColorBrewer",
  "shiny",
  "plotly",
  "sortable"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
    message("Installed: ", p)
  } else {
    message("Already installed: ", p)
  }
}

message("All packages ready.")
