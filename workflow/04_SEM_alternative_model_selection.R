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

data_rec <- data %>% filter(Treatment3 != "old-growth forest")

vars_to_scale <- c("ConIndex", "StrIndex", "multi_pol", "multi_sd", "FDSeeds", "FDSdlng")

data_scaled <- data_rec
data_scaled[vars_to_scale] <- lapply(data_rec[vars_to_scale], scale)

complete_rows <- complete.cases(data_scaled[,c("multi_pol", "multi_sd", "ConIndex", "StrIndex", "FDSeeds", "FDSdlng")]) 
complete_data <- data_scaled[complete_rows,c("multi_pol", "multi_sd", "ConIndex", "StrIndex", "FDSeeds", "FDSdlng")]

bf_pol_A   <- bf(multi_pol ~ ConIndex + StrIndex)
bf_sd_A    <- bf(multi_sd ~ ConIndex + StrIndex)
bf_seeds_A <- bf(FDSeeds | mi() ~ multi_pol + multi_sd + ConIndex + StrIndex)
bf_sdlng_A <- bf(FDSdlng | mi() ~ mi(FDSeeds) + ConIndex + StrIndex)

fit_full <- brm(
  bf_pol_A + bf_sd_A + bf_seeds_A + bf_sdlng_A + set_rescor(FALSE),
  data = data_scaled, cores = 4, chains = 4, iter = 4000
)

summary(fit_full)

# Supondo que você tenha 8 preditores totais no sistema e ache que 3 são reais
p_total <- 11
p_esperado <- 3
ratio <- p_esperado / (p_total - p_esperado)

# Definindo o prior
# df = 1 é o padrão, mas para N pequeno, df_global = 1 ajuda na estabilidade
# Rode isso para ver a lista de nomes que o brms aceita para priors
get_prior(bf_pol_A + bf_sd_A + bf_seeds_A + bf_sdlng_A + set_rescor(FALSE), 
          data = data_scaled)
meu_prior <- 
  set_prior(horseshoe(df = 1, par_ratio = ratio), class = "b", resp = "multipol") +
  set_prior(horseshoe(df = 1, par_ratio = ratio), class = "b", resp = "multisd") +
  set_prior(horseshoe(df = 1, par_ratio = ratio), class = "b", resp = "FDSeeds") +
  set_prior(horseshoe(df = 1, par_ratio = ratio), class = "b", resp = "FDSdlng")
# Rodando o modelo
fit_regularized <- brm(
  bf_pol_A + bf_sd_A + bf_seeds_A + bf_sdlng_A + set_rescor(FALSE),
  data = data_scaled,
  prior = meu_prior,
  cores = 4, chains = 4, iter = 4000, # Aumentei o iter para garantir convergência
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  file = "output/SEM_models/fit_regularized"
)
summary(fit_regularized)
plot(fit_regularized)

m_pol <- brm(
  bf_pol_A,
  data = data_scaled,
  prior = set_prior(horseshoe(df = 1, par_ratio = ratio), class = "b"),
  cores = 4, chains = 4, iter = 4000, # Aumentei o iter para garantir convergência
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  file = "output/SEM_models/fit_pol_3_rec_reg"
)
ref_pol <- get_refmodel(m_pol, resp = "multipol")
cv_pol <- cv_varsel(ref_pol)
suggest_size(cv_pol)
ranking(cv_pol)

m_sd <- brm(
  bf_sd_A,
  data = data_scaled,
  prior = set_prior(horseshoe(df = 1, par_ratio = ratio), class = "b"),
  cores = 4, chains = 4, iter = 4000, # Aumentei o iter para garantir convergência
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  file = "output/SEM_models/fit_sd_3_rec_reg"
)
ref_sd <- get_refmodel(m_sd, resp = "multisd")
cv_sd <- cv_varsel(ref_sd)
suggest_size(cv_sd)
ranking(cv_sd)

m_seeds <- brm(
  bf(multi_pol ~ 1) +
    bf(multi_sd ~ StrIndex) +  
    bf_seeds_A,
  data = complete_data,
  prior = set_prior(horseshoe(df = 1, par_ratio = ratio), class = "b", resp = "multisd") +
    set_prior(horseshoe(df = 1, par_ratio = ratio), class = "b", resp = "FDSeeds"),
  cores = 4, chains = 4, iter = 4000, # Aumentei o iter para garantir convergência
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  file = "output/SEM_models/fit_seeds_15_rec_reg"
)
ref_seeds <- get_refmodel(m_seeds, resp = "FDSeeds")
cv_seeds <- cv_varsel(ref_seeds)
suggest_size(cv_seeds)
ranking(cv_seeds)
