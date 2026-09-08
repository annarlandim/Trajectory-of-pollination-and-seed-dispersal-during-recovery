
# Reproducible entry point for the project
# Usage: source("main_script.R")

# if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::activate()

# Helpful project paths
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

pkgs <- c(
  "dplyr",
  "tidyr",
  "stringr",
  "psych",
  "ggplot2",
  "patchwork",
  "purrr",
  "nimble",
  "coda",
  "lattice",
  "MCMCvis",
  "vegan", 
  "parallel", 
  "bootnet",
  "boot",
  "qgraph",
  "brms", # tb nao vou usar
  "projpred", # n vou usar
  "huge",
  "rjags",
  "flextable",
  "officer"
)

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)

lapply(pkgs, library, character.only = TRUE)

renv::snapshot()

# Source utility functions
utils_files <- list.files(here("R"), pattern = "\\.R$", full.names = TRUE)
lapply(utils_files, source)

# Run workflow steps
source(here("workflow", "01_data_preparation.R"))
source(here("workflow", "02_data_analysis.R"))

message("Workflow completed.")

