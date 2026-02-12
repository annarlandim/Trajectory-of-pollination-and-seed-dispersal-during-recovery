### BRMS SEM ####

data <- read.csv("data/processed/model_df.csv")

weights_df <- read.csv("data/processed/weights_df.csv")

w_pol <- c(FDBees=weights_df[1,"prctg"]/sum(weights_df[1:3,"prctg"]), FDMoths=weights_df[2,"prctg"]/sum(weights_df[1:3,"prctg"]), FDBat_pol=weights_df[3,"prctg"]/sum(weights_df[1:3,"prctg"])) 
w_sd <- c(FDBats=weights_df[4,"prctg"]/sum(weights_df[4:6,"prctg"]), FDBirds=weights_df[5,"prctg"]/sum(weights_df[4:6,"prctg"]), FDNf=weights_df[6,"prctg"]/sum(weights_df[4:6,"prctg"]))

data$multi_pol <- apply(data[, names(w_pol)], 1, function(x) {
  sum(x * w_pol, na.rm = TRUE)
})
data$multi_sd <-apply(data[, names(w_sd)], 1, function(x) {
  sum(x * w_sd, na.rm = TRUE)
})

data_rec <- data %>% filter(Treatment3 != "old-growth forest")

vars_to_scale <- c("ConIndex", "StrIndex", "multi_pol", "multi_sd", "FDSeeds", "FDSdlng")

data_scaled <- data
data_scaled[vars_to_scale] <- lapply(data[vars_to_scale], scale)

# --- Submodels ---

bf_pol <- bf(
  multi_pol ~ ConIndex + StrIndex
)

bf_sd <- bf(
  multi_sd ~ ConIndex + StrIndex
)

bf_seeds <- bf(FDSeeds | mi() ~ ConIndex + StrIndex + multi_pol + multi_sd)

bf_sdlng <- bf(FDSdlng | mi() ~ StrIndex + mi(FDSeeds))

# --- Combine as SEM ---
bayes_sem <- brm(
  bf_pol + bf_sd + bf_seeds + bf_sdlng +
    set_rescor(FALSE),   # matches your DAG: no unexplained covariances
  data = data_scaled,
  chains = 4,
  iter = 4000,
  cores = 4,
  control = list(adapt_delta = 0.95)
)

summary(bayes_sem)

# Verifique visualmente os caminhos
conditional_effects(bayes_sem, "StrIndex", resp = "FDSdlng")
conditional_effects(bayes_sem, "FDSeeds", resp = "FDSdlng")

## Mediated effects ##

draws <- as_draws_df(bayes_sem) %>%
  mutate(
    
    # ----- MEDIATION TO SEEDS -----
    
    ind_Str_multisd_Seeds = b_multisd_StrIndex * b_FDSeeds_multi_sd,
    
    ind_Str_multipol_Seeds = b_multipol_StrIndex * b_FDSeeds_multi_pol,
    
    total_Str_Seeds =
      b_FDSeeds_StrIndex +
      ind_Str_multisd_Seeds +
      ind_Str_multipol_Seeds,
    
    
    # ----- MEDIATION TO SEEDLINGS -----
    
    ind_Str_Seeds_Sdlng = b_FDSeeds_StrIndex * bsp_FDSdlng_miFDSeeds,
    
    total_Str_Sdlng =
      b_FDSdlng_StrIndex +
      ind_Str_Seeds_Sdlng,
    
    ind_Pol_Seeds_Sdlng =
      b_FDSeeds_multi_pol * bsp_FDSdlng_miFDSeeds,
    
    ind_SD_Seeds_Sdlng =
      b_FDSeeds_multi_sd * bsp_FDSdlng_miFDSeeds
  )

draws %>%
  pivot_longer(
    cols = starts_with(c("ind_", "total_")),
    names_to = "effect",
    values_to = "value"
  ) %>%
  group_by(effect) %>%
  summarise(
    mean = mean(value),
    lower = quantile(value, .025),
    upper = quantile(value, .975),
    prob_positive = mean(value > 0)
  )
