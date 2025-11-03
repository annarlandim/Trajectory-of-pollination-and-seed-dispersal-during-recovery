### Organizing project

# How to use daily
# 
# Do your edits in R/, workflow/, and save outputs to output/.
# 
# Run everything via:
#   
#   source("main_script.R")
# 
# 
# When you add/remove packages:
#   
#   renv::snapshot()
# git add renv.lock
# git commit -m "Update deps"

#######
dirs <- c(
  "data/raw", "data/processed",
  "R",                 # your utility functions
  "workflow",          # scripts called by main_script
  "output",
  "renv"               # will be populated by renv
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

write_if_absent <- function(path, txt) if (!file.exists(path)) writeLines(txt, path)

write_if_absent("main_script.R", '
# Reproducible entry point for the project
# Usage: source("main_script.R")

# Load/activate renv (will auto-install on first run if needed)
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::activate()

# Helpful project paths
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# Source utility functions
utils_files <- list.files(here("R"), pattern = "\\\\.R$", full.names = TRUE)
invisible(lapply(utils_files, source))

# Run workflow steps
source(here("workflow", "01_data_preparation.R"))
source(here("workflow", "02_data_analysis.R"))

message("Workflow completed.")
')

rite_if_absent("workflow/01_data_preparation.R", '
# Read raw data, clean, and write processed data
# Example:
# raw <- read.csv(here::here("data", "raw", "your_raw_file.csv"))
# processed <- raw  # <- do your cleaning here
# saveRDS(processed, here::here("data", "processed", "processed.rds"))
')

# analysis stub
write_if_absent("workflow/02_data_analysis.R", '
# Read processed data and produce results/figures/tables
# processed <- readRDS(here::here("data", "processed", "processed.rds"))
# ... analysis ...
# ggplot2, models, outputs saved to output/
')

# utils stub
write_if_absent("R/utils.R", '
# Put helper functions here, they are sourced by main_script.R
# example:
# not_na <- function(x) x[!is.na(x)]
')

# README
write_if_absent("README.md", "# Multifunctionality_recovery_time\n\nReproducible R project scaffold with renv.\n\n## How to run\n1. Open the .Rproj in RStudio\n2. Run: `source(\"main_script.R\")`\n")

# # .gitignore (keeps lockfile tracked, ignores renv library binaries)
# write_if_absent(".gitignore", '
# .Rhistory
# .RData
# .Ruserdata
# .Rproj.user/
# .Rproj.user
# .Rproj.user/*
# .Rproj.user/*/shared/notebooks/*
# .Rproj.user/*/suspended-data
# .Rproj.user/*/pcs/*
# .Rproj.user/*/sources/*
# .Rproj.user/*/sessions/*
# .Rproj.user/*/profiles/*
# .Rproj.user/*/explorer-cache/*
# .Rproj.user/shared/notebooks/paths
# .Rproj.user/shared/notebooks/indices
# .Rproj.user/shared/notebooks/lock_file
# .Rproj.user/shared/notebooks/*
# .Rproj.user/sources/prop/*
# 
# # renv cache is global; keep lockfile tracked but ignore local library
# renv/library/
# renv/staging/
# 
# # data that shouldn’t be versioned (optional)
# data/processed/
# output/
# ')

# Optional: RStudio project options — newline at end of file etc. (safe default)
# write_if_absent("Multifunctionality_recovery_time.Rproj", '
# Version: 1.0
# 
# RestoreWorkspace: No
# SaveWorkspace: No
# AlwaysSaveHistory: Default
# 
# EnableCodeIndexing: Yes
# UseSpacesForTab: Yes
# NumSpacesForTab: 2
# Encoding: UTF-8
# 
# RnwWeave: knitr
# LaTeX: pdfLaTeX
# ')

#######

install.packages("renv")    # first time only
renv::init()                # creates renv.lock and renv/infra
