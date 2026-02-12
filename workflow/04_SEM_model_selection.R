set.seed(1212)

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

# data_rec <- data %>% filter(Treatment3 != "old-growth forest")

vars_to_scale <- c("ConIndex", "StrIndex", "multi_pol", "multi_sd", "FDSeeds", "FDSdlng")

data_scaled <- data
data_scaled[vars_to_scale] <- lapply(data[vars_to_scale], scale)

## Models:

# Saturated (A)
# bf_pol_A   <- bf(multi_pol ~ ConIndex + StrIndex)
# bf_sd_A    <- bf(multi_sd ~ ConIndex + StrIndex)
# mi() allows missing data in seeds and seedlings:
# bf_seeds_A <- bf(FDSeeds | mi() ~ ConIndex + StrIndex + multi_pol + multi_sd)
# bf_sdlng_A <- bf(FDSdlng | mi() ~ ConIndex + StrIndex + mi(FDSeeds)) 
# 
# fit_saturated <- brm(
#   bf_pol_A + bf_sd_A + bf_seeds_A + bf_sdlng_A + set_rescor(FALSE),
#   data = data_scaled, cores = 4, chains = 4, iter = 4000,
#   file = "output/SEM_models/fit_saturated" # Salva o modelo para não precisar rodar sempre
# )

fit_saturated <- readRDS("output/SEM_models/fit_saturated.RDS")

summary(fit_saturated)

# Mediated: Str and con only via pol and sd (B)
# bf_pol_B   <- bf(multi_pol ~ ConIndex + StrIndex)
# bf_sd_B    <- bf(multi_sd ~ ConIndex + StrIndex)
# bf_seeds_B <- bf(FDSeeds | mi() ~ multi_pol + multi_sd) # Sem Con e Str diretos
# bf_sdlng_B <- bf(FDSdlng | mi() ~ mi(FDSeeds))         # Só sementes importam para plântulas
# 
# fit_mediated <- brm(
#   bf_pol_B + bf_sd_B + bf_seeds_B + bf_sdlng_B + set_rescor(FALSE),
#   data = data_scaled, cores = 4, chains = 4, iter = 4000,
#   file = "output/SEM_models/fit_mediated"
# )

fit_mediated <- readRDS("output/SEM_models/fit_mediated.RDS")

summary(fit_mediated)

# Ab+Bio: Str and biotic affect seedlings (C)
# bf_pol_C   <- bf(multi_pol ~ ConIndex + StrIndex)
# bf_sd_C    <- bf(multi_sd ~ ConIndex + StrIndex)
# bf_seeds_C <- bf(FDSeeds | mi() ~ multi_pol + multi_sd) 
# bf_sdlng_C <- bf(FDSdlng | mi() ~ mi(FDSeeds) + StrIndex)         
# 
# fit_abbio <- brm(
#   bf_pol_C + bf_sd_C + bf_seeds_C + bf_sdlng_C + set_rescor(FALSE),
#   data = data_scaled, cores = 4, chains = 4, iter = 4000,
#   file = "output/SEM_models/fit_abbio"
# )

fit_abbio <- readRDS("output/SEM_models/fit_abbio.RDS")

summary(fit_abbio)

# Ab: Str affect seedlings (E)
# bf_pol_E   <- bf(multi_pol ~ ConIndex + StrIndex)
# bf_sd_E    <- bf(multi_sd ~ ConIndex + StrIndex)
# bf_seeds_E <- bf(FDSeeds | mi() ~ multi_pol + multi_sd)
# bf_sdlng_E <- bf(FDSdlng | mi() ~ StrIndex)
# 
# fit_ab <- brm(
#   bf_pol_E + bf_sd_E + bf_seeds_E + bf_sdlng_E + set_rescor(FALSE),
#   data = data_scaled, cores = 4, chains = 4, iter = 4000,
#   file = "output/SEM_models/fit_ab"
# )

fit_ab <- readRDS("output/SEM_models/fit_ab.RDS")

summary(fit_ab)

# Null: Nothing affects seedlings (D)
# bf_pol_D   <- bf(multi_pol ~ ConIndex + StrIndex)
# bf_sd_D    <- bf(multi_sd ~ ConIndex + StrIndex)
# bf_seeds_D <- bf(FDSeeds | mi() ~ multi_pol + multi_sd) 
# bf_sdlng_D <- bf(FDSdlng | mi() ~ 1)         
# 
# fit_null <- brm(
#   bf_pol_D + bf_sd_D + bf_seeds_D + bf_sdlng_D + set_rescor(FALSE),
#   data = data_scaled, cores = 4, chains = 4, iter = 4000,
#   file = "output/SEM_models/fit_null"
# )

fit_null <- readRDS("output/SEM_models/fit_null.RDS")

summary(fit_null)

# to allow NAs in model comparison:
# 
# # brms divides the 62 plots in 10 random groups 
# kfold_sat <- kfold(fit_saturated, K = 10, resp = c("multipol", "multisd"), chains = 4, seed = 1212) 
# kfold_med <- kfold(fit_mediated, K = 10, resp = c("multipol", "multisd"), chains = 4, seed = 1212)  
# kfold_abbio <- kfold(fit_abbio, K = 10, resp = c("multipol", "multisd"), chains = 4, seed = 1212) 
# kfold_null <- kfold(fit_null, K = 10, resp = c("multipol", "multisd"), chains = 4, seed = 1212) 

# check if the absolute difference |(elpd_diff)| is two times larger than se (se_diff).
# if not, models are the same
# comp_kfold <- loo_compare(kfold_sat, kfold_med, kfold_abbio, kfold_null)
# comp_kfold

complete_rows <- complete.cases(data_scaled[,c("multi_pol", "multi_sd", "ConIndex", "StrIndex", "FDSeeds", "FDSdlng")]) 
complete_data <- data_scaled[complete_rows,c("multi_pol", "multi_sd", "ConIndex", "StrIndex", "FDSeeds", "FDSdlng")]

loo_sat <- loo(fit_saturated, newdata = complete_data)
loo_med <- loo(fit_mediated, newdata = complete_data)
loo_abbio <- loo(fit_abbio, newdata = complete_data)
loo_ab <- loo(fit_ab, newdata = complete_data)
loo_null <- loo(fit_null, newdata = complete_data)

comp_loo <- loo_compare(loo_sat, loo_med, loo_abbio, loo_ab, loo_null)

## Final diagnostics:

pp_check(fit_ab, resp = "multipol")
pp_check(fit_ab, resp = "multisd")
pp_check(fit_ab, resp = "FDSeeds")
pp_check(fit_ab, resp = "FDSdlng")

bayes_R2(fit_ab, na.rm = T)
conditional_effects(fit_ab, resp = "FDSdlng")

# Função para calcular R2 de uma variável específica em modelos multivariados
get_r2_manual <- function(model, resp_name, sigma_name) {
  # Extrai as predições (linear predictor)
  mu <- fitted(model, resp = resp_name, summary = FALSE)
  var_mu <- apply(mu, 1, var)
  
  # Extrai o sigma do posterior
  sigma <- as_draws_df(model)[[sigma_name]]
  var_epsilon <- sigma^2
  
  # Calcula o R2 para cada draw
  r2_draws <- var_mu / (var_mu + var_epsilon)
  
  # Retorna o resumo (média e intervalos)
  return(data.frame(
    Estimate = mean(r2_draws),
    Q2.5 = quantile(r2_draws, 0.025),
    Q97.5 = quantile(r2_draws, 0.975)
  ))
}

# Exemplos de uso:
r2_seeds <- get_r2_manual(fit_abbio, "FDSeeds", "sigma_FDSeeds")
r2_sdlng <- get_r2_manual(fit_abbio, "FDSdlng", "sigma_FDSdlng")

# checking if it works:
r2_pol <- get_r2_manual(fit_abbio, "multipol", "sigma_multipol")
r2_sd <- get_r2_manual(fit_abbio, "multisd", "sigma_multisd")
