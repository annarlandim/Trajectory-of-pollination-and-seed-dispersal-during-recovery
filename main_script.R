
# Reproducible entry point for the project
# Usage: source("main_script.R")

if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
# renv::activate()
renv::restore()

library(here)

# Source utility functions
utils_files <- list.files(here("R"), pattern = "\\.R$", full.names = TRUE)
lapply(utils_files, source)

# Run workflow steps
source(here("workflow", "01_Data_cleaning.R"))
source(here("workflow", "02_Model_data_prep.R"))
source(here("workflow", "03_Bayesian_models.R"))

message("Workflow completed.")

