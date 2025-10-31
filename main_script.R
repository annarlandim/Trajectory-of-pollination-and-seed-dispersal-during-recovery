
# Reproducible entry point for the project
# Usage: source("main_script.R")

# Load/activate renv (will auto-install on first run if needed)
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::activate()

# Helpful project paths
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# Source utility functions
utils_files <- list.files(here("R"), pattern = "\\.R$", full.names = TRUE)
invisible(lapply(utils_files, source))

# Run workflow steps
source(here("workflow", "01_data_preparation.R"))
source(here("workflow", "02_data_analysis.R"))

message("Workflow completed.")

