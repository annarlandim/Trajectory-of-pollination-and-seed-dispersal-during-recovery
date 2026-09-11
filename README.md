#################
DATA AND CODE FROM
#################

"Divergent trajectories in the functional diversity of pollination and seed-dispersal interactions during forest recovery"

#######
CONTENT
#######

## R scripts — main pipeline (run in order via main_script.R)

1) main_script.R - Reproducible entry point. Restores the R package environment via `renv` and sources the three workflow scripts below in order.
2) workflow/01_Data_cleaning.R - Cleans and harmonizes raw interaction, trait, and connectivity data; standardizes plant and animal species names; writes processed interaction and trait tables.
3) workflow/02_Model_data_prep.R - Builds trait spaces (PCA) per function, computes per-plot functional diversity metrics per group, and assembles the final model data frame.
4) workflow/03_Bayesian_models.R - Fits the hierarchical NIMBLE model estimating  trajectories and the effect of forest connectivity.

## R scripts — figures, diagnostics, and checks (run separately)

5) workflow/02_01_Figure_trait_spaces.R - Produces the trait-space figures (overall and per recovery stage) shown in the manuscript.
6) workflow/03_00_Priors_sanity_check.R - Prior predictive checks for the priors used in the Bayesian model.
7) workflow/03_01_Figures_trajectories_and_net_change.R - Produces the trajectory and net-change figures.
8) workflow/03_02_Model_diagnostics.R - Summarizes MCMC diagnostics (effective sample size, PSRF, credible intervals) into a table.

## Raw data (data/raw/)

Interaction, trait, and connectvity data collected in the field and compiled from prior REASSEMBLY project datasets (subprojects SP3: pollination; SP4: seed dispersal; and Central Management), including:
- Plot-level treatment and regeneration-time metadata (Central Management)
- Pollination interaction and trait records (SP3/)
- Seed-dispersal interaction and trait records, including camera-trap, direct-observation, and mist-net (bat) data (SP4/)
- Connectivity data (Central Management)

## Processed data (data/processed/, produced by the scripts)

Cleaned interaction tables and trait tables per interaction group (bees, moths, bats-pollination, bats-seed dispersal, birds, non-flying mammals), the assembled connectivity index, and the final model data frame used as input to the Bayesian models.

## Model outputs (output/, produced by the scripts)

Posterior samples,  model metadata, and diagnostic summaries from the Bayesian models, and the SVG figures reproduced in the manuscript.

############
REQUIREMENTS
############

- R (version: 4.5.3 (2026-03-11 ucrt))
- `renv` for package version management — run `renv::restore()` before anything else
- A working C++ compiler toolchain for NIMBLE to compile models: Rtools (Windows) or Xcode Command Line Tools (Mac)

############
SESSION INFO
############

R version 4.5.3 (2026-03-11 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=English_Germany.utf8  LC_CTYPE=English_Germany.utf8    LC_MONETARY=English_Germany.utf8
[4] LC_NUMERIC=C                     LC_TIME=English_Germany.utf8    

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.5.3    tools_4.5.3       rstudioapi_0.18.0 renv_1.1.5