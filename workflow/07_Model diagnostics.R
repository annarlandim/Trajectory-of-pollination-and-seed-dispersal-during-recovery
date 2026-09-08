#### MCMC diagnostics summary ####
# Consolidates Gelman-Rubin (R-hat) and effective sample size (ESS)
# across all monitored parameter blocks into one table, for reporting
# in the Methods / supplementary material.

library(coda)
library(dplyr)
library(purrr)

# --- Load your saved posterior samples ---
# Adjust path to wherever you saved the mcmc.list object
# (the commented-out line in your original script:
#  saveRDS(samples, "output/model_posteriors.4_con.rds"))

samples <- readRDS("output/model_posteriors_con.rds")
samples_all <- as.mcmc(do.call(rbind, samples))

# --- Parameter blocks to check ---
# 'log_scale = TRUE' for parameters that are strictly positive and were
# modelled on the log scale (theta_inf, theta_0), matching your existing
# gelman_log() approach. Hyperparameters and connectivity slopes/intercepts
# are already unconstrained, so they're checked on their native scale.
#
# 'scope' controls which numeric indices are kept, since this saved object
# includes seeds and seedlings groups that are NOT used in the paper:
#   - "taxon": per-interaction-group parameters (theta_inf, theta_0,
#     alpha_con, beta_con). Indices 1:3 = pollination groups (bees, moths,
#     bat_pol), 4:6 = seed-dispersal groups (bats, birds, nf), 7:10 = seeds,
#     11:14 = seedlings. Only 1:6 are kept.
#   - "hypergroup": function-level hyperparameters (mean_*/sigma_* blocks).
#     Index 1 = pollination, 2 = seed dispersal, 3 = seeds, 4 = seedlings.
#     Only 1:2 are kept.

param_blocks <- list(
  list(name = "theta_inf",                  pattern = "^theta_inf\\[",                         log_scale = TRUE,  scope = "taxon"),
  list(name = "theta_0",                    pattern = "^theta_0\\[",                           log_scale = TRUE,  scope = "taxon"),
  list(name = "alpha_con",                  pattern = "^alpha_con\\[",                         log_scale = FALSE, scope = "taxon"),
  list(name = "beta_con",                   pattern = "^beta_con\\[",                          log_scale = FALSE, scope = "taxon"),
  list(name = "mean_theta_inf/sigma_theta_inf", pattern = "^(mean_theta_inf|sigma_theta_inf)\\[", log_scale = FALSE, scope = "hypergroup"),
  list(name = "mean_theta_0/sigma_theta_0",     pattern = "^(mean_theta_0|sigma_theta_0)\\[",     log_scale = FALSE, scope = "hypergroup"),
  list(name = "mean_alpha_con/sigma_alpha_con", pattern = "^(mean_alpha_con|sigma_alpha_con)\\[", log_scale = FALSE, scope = "hypergroup"),
  list(name = "mean_beta_con/sigma_beta_con",   pattern = "^(mean_beta_con|sigma_beta_con)\\[",   log_scale = FALSE, scope = "hypergroup")
)

# IMPORTANT: this assumes the k-index ordering (1 = pollination, 2 = seed
# dispersal, 3 = seeds, 4 = seedlings) is the same across ALL hypergroup
# blocks, since they should all derive from the same group_index_vec.
# Verify this holds for mean_theta_inf/sigma_theta_inf,
# mean_alpha_con/sigma_alpha_con, and mean_beta_con/sigma_beta_con -
# not just mean_theta_0/sigma_theta_0 - before trusting the filtered output.

# --- Helper: get R-hat + ESS for one block ---

diagnose_block <- function(block, samples, samples_all) {
  
  names_vars <- varnames(samples)
  idx <- grep(block$pattern, names_vars)
  
  if (length(idx) == 0) {
    warning(sprintf("No parameters matched for block '%s'", block$name))
    return(NULL)
  }
  
  if (block$log_scale) {
    # Gelman on log scale (as in your gelman_log() function)
    samples_list <- as.mcmc.list(lapply(samples, function(chain) {
      as.mcmc(log(chain[, idx, drop = FALSE]))
    }))
    ess_input <- log(samples_all[, idx, drop = FALSE])
  } else {
    samples_list <- samples[, idx, drop = FALSE]
    ess_input <- samples_all[, idx, drop = FALSE]
  }
  
  gd <- gelman.diag(samples_list, autoburnin = FALSE, multivariate = FALSE)
  rhat_vals <- gd$psrf[, "Point est."]
  
  ess_vals <- effectiveSize(ess_input)
  
  tibble(
    block = block$name,
    parameter = names(rhat_vals),
    rhat = as.numeric(rhat_vals),
    ess = as.numeric(ess_vals[names(rhat_vals)])
  )
}

# --- Run across all blocks ---

diag_table_all <- map_dfr(param_blocks, diagnose_block, samples = samples, samples_all = samples_all)

# --- Filter out seeds/seedlings parameters, not used in this paper ---

diag_table <- diag_table_all %>%
  mutate(
    param_index = as.integer(str_extract(parameter, "(?<=\\[)\\d+(?=\\])")),
    scope = map_chr(block, ~ param_blocks[[which(map_chr(param_blocks, "name") == .x)]]$scope)
  ) %>%
  filter(
    (scope == "taxon"      & param_index <= 6) |
      (scope == "hypergroup" & param_index <= 2)
  ) %>%
  select(-scope)

# --- Full table (per-parameter, pollination + seed dispersal only) ---
print(diag_table, n = Inf)

# For full transparency / a sanity check, this is what was excluded:
excluded <- anti_join(diag_table_all, diag_table, by = c("block", "parameter"))
cat(sprintf("Excluded %d seeds/seedlings parameters from reporting.\n", nrow(excluded)))

# --- Summary for reporting in text ---
# This gives you the numbers to fill into:
# "R-hat < [X] for all parameters, ESS > [X]"

summary_stats <- diag_table %>%
  summarise(
    n_parameters = n(),
    max_rhat = max(rhat, na.rm = TRUE),
    min_ess  = min(ess,  na.rm = TRUE),
    n_rhat_above_1.01 = sum(rhat > 1.01, na.rm = TRUE),
    n_rhat_above_1.05 = sum(rhat > 1.05, na.rm = TRUE),
    n_ess_below_400   = sum(ess < 400, na.rm = TRUE)
  )

print(summary_stats)

# --- Optional: flag any problem parameters individually ---
# Useful if max_rhat or min_ess look concerning and you need to find
# which specific group/parameter is driving it, rather than assuming
# it's uniform across all six interaction groups.

problem_params <- diag_table %>%
  filter(rhat > 1.01 | ess < 400) %>%
  arrange(desc(rhat))

print(problem_params)

#### Supplementary diagnostics table ####
# Builds a publication-ready table of posterior summaries and convergence
# diagnostics (PSRF, Neff) for each interaction group, matching the style
# of the reference table. Exports directly to a .docx file.

library(dplyr)
library(purrr)
library(stringr)
library(coda)
library(flextable)
library(officer)

samples <- readRDS("output/model_posteriors_con.rds")
samples_all <- as.mcmc(do.call(rbind, samples))

# --- Interaction group labels ---
# Order matches column order in dataSub: FDBees, FDMoths, FDBat_pol,
# FDBats, FDBirds, FDNf. Adjust labels if your naming differs.

group_labels <- c(
  "1" = "Bees",
  "2" = "Moths",
  "3" = "Bats (pollination)",
  "4" = "Bats (seed dispersal)",
  "5" = "Birds",
  "6" = "Non-flying mammals"
)

# --- Helper: quantiles + PSRF + Neff for one parameter type across groups ---
# transform_fn is applied to each posterior draw before summarizing, so it's
# used to convert the monitored sigma_old/sigma_rec into precision (tau).

summarise_param <- function(base_name, pretty_name, samples, samples_all,
                            indices = 1:6, transform_fn = identity) {
  
  pattern <- paste0("^", base_name, "\\[")
  names_vars <- varnames(samples)
  idx_cols <- grep(pattern, names_vars)
  
  col_idx_num <- as.integer(str_extract(names_vars[idx_cols], "(?<=\\[)\\d+(?=\\])"))
  keep <- idx_cols[col_idx_num %in% indices]
  keep <- keep[order(col_idx_num[col_idx_num %in% indices])]
  
  draws <- samples_all[, keep, drop = FALSE]
  draws_t <- apply(draws, 2, transform_fn)
  
  q <- t(apply(draws_t, 2, quantile, probs = c(0.5, 0.025, 0.25, 0.75, 0.975)))
  colnames(q) <- c("Median", "q2.5", "q25", "q75", "q97.5")
  
  neff <- effectiveSize(as.mcmc(draws_t))
  
  samples_t_list <- as.mcmc.list(lapply(samples, function(chain) {
    sub <- chain[, keep, drop = FALSE]
    as.mcmc(apply(sub, 2, transform_fn))
  }))
  psrf <- gelman.diag(samples_t_list, autoburnin = FALSE, multivariate = FALSE)$psrf[, "Point est."]
  
  tibble(
    Parameter = pretty_name,
    Group = group_labels[as.character(sort(indices))],
    Median = q[, "Median"], `2.5%` = q[, "q2.5"], `25%` = q[, "q25"],
    `75%` = q[, "q75"], `97.5%` = q[, "q97.5"],
    Neff = as.numeric(neff), PSRF = as.numeric(psrf)
  )
}

tau <- function(sigma) 1 / sigma^2

# --- Build table: one block per parameter, six groups each ---
# Greek/subscript labels match the reference table (Theta0, ThetaOG, alpha,
# beta, TauOG, TauREC). theta_inf = the old-growth asymptote (ThetaOG);
# theta_0 = the active-land-use starting value (Theta0).

table_long <- bind_rows(
  summarise_param("theta_0",   "\u03980",   samples, samples_all),
  summarise_param("theta_inf", "\u0398OG",  samples, samples_all),
  summarise_param("alpha_con", "\u03b1",    samples, samples_all),
  summarise_param("beta_con",  "\u03b2",    samples, samples_all),
  summarise_param("sigma_old", "\u03a4OG",  samples, samples_all, transform_fn = tau),
  summarise_param("sigma_rec", "\u03a4REC", samples, samples_all, transform_fn = tau)
)

# --- Order rows: grouped by interaction group, cycling through parameters ---
param_order <- c("\u03980", "\u0398OG", "\u03b1", "\u03b2", "\u03a4OG", "\u03a4REC")
group_order <- group_labels

table_final <- table_long %>%
  mutate(
    Parameter = factor(Parameter, levels = param_order),
    Group = factor(Group, levels = group_order)
  ) %>%
  arrange(Group, Parameter) %>%
  mutate(
    across(c(Median, `2.5%`, `25%`, `75%`, `97.5%`, PSRF), ~ round(.x, 2)),
    Neff = round(Neff, 1),
    `Parameter ` = paste(Parameter, Group)   # combined label, matching reference table style
  ) %>%
  select(`Parameter ` , Median, `2.5%`, `25%`, `75%`, `97.5%`, Neff, PSRF)

# --- Build Word-ready table ---

ft <- flextable(table_final) %>%
  set_header_labels(`Parameter ` = "Parameter") %>%
  autofit() %>%
  theme_booktabs() %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  add_footer_lines(
    "PSRF: potential scale reduction factor (values close to 1 indicate convergence). Neff: effective sample size. \u03980: functional diversity in active land use. \u0398OG: functional diversity in old-growth forest. \u03b1, \u03b2: intercept and slope of the recovery rate. \u03a4OG, \u03a4REC: precision of observations in old-growth and recovering forests, respectively."
  ) %>%
  fontsize(size = 8, part = "footer")

# --- Export to Word ---
save_as_docx(ft, path = "output/TableS1_diagnostics.docx")

# Also save as CSV for your own reference / re-use
write.csv(table_final, "output/TableS1_diagnostics.csv", row.names = FALSE)
