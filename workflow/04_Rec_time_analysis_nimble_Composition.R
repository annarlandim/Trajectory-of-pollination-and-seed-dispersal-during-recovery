library(dplyr)
library(tidyr)
library(purrr)

library(nimble)
library(coda)
library(lattice)
library(MCMCvis)


#### Data ####

data <- read.csv("data/processed/model_df.csv")

weights_df <- read.csv("data/processed/weights_df.csv")

data <- data %>% mutate(ConIndex = scale(ConIndex))

data$type <- factor(ifelse(1:nrow(data) %in% grep("OG", data$Plot_ID), "old", "rec"),
                    levels = c("old", "rec"))

dataSub <- subset(data, select = c(type, RegTime, ConIndex,  
                                   FDBees, FC1Bees, FC2Bees, FDMoths, FC1Moths, FC2Moths, FDBat_pol, FC1Bat_pol, FC2Bat_pol,
                                   FDBats, FC1Bats, FC2Bats, FDBirds, FC1Birds, FC2Birds, FDNf, FC1Nf, FC2Nf,
                                   VerticalVH, MaxTH, AGB))

table(dataSub$type)
str(dataSub)


# plot the data

par(mfrow = c(1, 3))
plot(FC1Bees ~ RegTime, dataSub)
plot(FC1Moths ~ RegTime, dataSub)
plot(FC1Bat_pol ~ RegTime, dataSub)
dev.off()

par(mfrow = c(1, 3))
plot(FC2Bees ~ RegTime, dataSub)
plot(FC2Moths ~ RegTime, dataSub)
plot(FC2Bat_pol ~ RegTime, dataSub)
dev.off()

par(mfrow = c(1, 3))
plot(FC1Bats ~ RegTime, dataSub)
plot(FC1Birds ~ RegTime, dataSub)
plot(FC1Nf ~ RegTime, dataSub)
dev.off()

par(mfrow = c(1, 3))
plot(FC2Bats ~ RegTime, dataSub)
plot(FC2Birds ~ RegTime, dataSub)
plot(FC2Nf ~ RegTime, dataSub)
dev.off()

# define data for jags model

# remove in the future. problem with (log(negative))
dataSub <- dataSub %>%
  select(-FDBees, -FDMoths, -FDBat_pol, -FDBats, -FDBirds, -FDNf, -VerticalVH, -MaxTH, -AGB) %>%
  mutate(across(
    .cols = -c(1:3),
    .fns = ~ as.numeric(scale(.x))
  ))

zeros <- dataSub %>% 
  map_lgl(~ any(. == 0, na.rm = TRUE))
zeros

long <- dataSub %>%
  pivot_longer(cols = -c(1:3), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = colnames(dataSub)[-c(1:3)]), variable = as.integer(variable)) %>%
  filter(!is.na(value))

sd_old_vec <- long %>%
  filter(type == "old") %>%
  group_by(variable) %>%
  summarise(sd = sd(value, na.rm = TRUE), .groups = "drop") %>%
  pull(sd)

# to avoid very tiny sds:
sd_old_vec <- pmax(sd_old_vec, 0.2)

sd_rec_vec <- long %>%
  filter(type == "rec") %>%
  group_by(variable) %>%
  summarise(sd = sd(value, na.rm = TRUE), .groups = "drop") %>%
  pull(sd)

old_df <- long %>%
  filter(type == "old") %>%
  transmute(Y_old = value,
            variable_old = variable) 

rec_df <- long %>%
  filter(type == "rec") %>%
  transmute(
    Y_rec = value,
    variable_rec = variable,
    tx = RegTime,
    connectivity = ConIndex
  )

group_index_vec <- c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2)

constList <- list(
  n_var = ncol(dataSub) - 3, # removing plot_id, reg_time and con_index
  n_groups = 2,          # mudar quando add as plantulas
  group_index = group_index_vec,
  n_old = nrow(old_df),
  n_rec = nrow(rec_df),
  variable_old = old_df$variable_old,
  variable_rec = rec_df$variable_rec
)

dataList <- list(
  Y_old = old_df$Y_old,
  s_old = sd_old_vec,
  Y_rec = rec_df$Y_rec,
  s_rec = sd_rec_vec,
  tx = rec_df$tx,
  connectivity = as.numeric(rec_df$connectivity)
)

#### Model ####

code_gaussian_stratified <- nimbleCode({
  
  ################
  #  Hyperpriors # 
  ################
  
  for(k in 1:n_groups){
    
    # --- theta_inf ---
    mean_theta_inf[k]  ~ dnorm(0, sd = 10)
    sigma_theta_inf[k] ~ dexp(1)
    
    # --- theta_0 ---
    mean_theta_0[k]    ~ dnorm(0, sd = 10)
    sigma_theta_0[k]   ~ dexp(1)
    
    # --- Connectivity (alpha and beta) ---
    mean_alpha_con[k]  ~ dnorm(0, sd = 1)
    sigma_alpha_con[k] ~ dexp(1)
    
    mean_beta_con[k]   ~ dnorm(0, sd = 1)
    sigma_beta_con[k]  ~ dexp(1)
  }
  
  ################################
  #   Parameters per variable j  #
  ################################
  
  # --- theta_inf ---
  
  for(j in 1:n_var){
    beta_inf_raw[j] ~ dnorm(0, sd = 1)
    # asymptote per group, original scale (always > 0)
    theta_inf[j] <- 
      mean_theta_inf[ group_index[j] ] + 
        beta_inf_raw[j] * sigma_theta_inf[ group_index[j] ]
    
    # priors on variance components:
    sigma_raw_old[j] ~ dexp(1)
    sigma_old[j] <- s_old[j] * sigma_raw_old[j] # scale by observed sd
    tau_old[j] <- 1/pow(sigma_old[j], 2)
  }
  
  # --- theta_0 ---
  
  for(j in 1:n_var){
    beta_0_raw[j] ~ dnorm(0, sd = 1)
    # asymptote per group, original scale (always > 0)
    theta_0[j] <- 
      mean_theta_0[ group_index[j] ] + 
        beta_0_raw[j] * sigma_theta_0[ group_index[j] ]
    
    # priors on variance components:
    sigma_raw_rec[j] ~ dexp(1)
    sigma_rec[j] <- s_rec[j] * sigma_raw_rec[j] # scale by observed sd
    tau_rec[j] <- 1/pow(sigma_rec[j], 2)
  }
  
  
  # --- alpha_con, beta_con ---
  
  for (j in 1:n_var) {
    alpha_con_raw[j] ~ dnorm(0, sd = 1)
    alpha_con[j] <- mean_alpha_con[ group_index[j] ] + 
      alpha_con_raw[j] * sigma_alpha_con[ group_index[j] ]
    
    beta_con_raw[j] ~ dnorm(0, sd = 1)
    beta_con[j] <- mean_beta_con[ group_index[j] ] + 
      beta_con_raw[j] * sigma_beta_con[ group_index[j] ]
  }
  
  ##########################
  # Likelihoods            #
  ##########################
  
  # old-growth
  for (i in 1:n_old) {
    # parameterization reminder:
    # meanlog = log(theta_inf[group]), sdlog = sigma_old[group]
    Y_old[i] ~ dnorm(mean = theta_inf[ variable_old[i] ],
                      sd = sigma_old[ variable_old[i] ])
  }
  
  # recovering forests
  for (i in 1:n_rec) {
    
    # recovery trajectories per group
    mu_rec[i] <-
      theta_0[ variable_rec[i] ] +
      (theta_inf[ variable_rec[i] ] - theta_0[ variable_rec[i] ]) *
      (1 - exp(-lambda[i] * tx[i]))
    
    # recovery rate (lambda) according to connectivity:
    lambda[i] <- exp(
      alpha_con[ variable_rec[i] ] +
        beta_con[  variable_rec[i] ] * connectivity[i]
    )
    
    Y_rec[i] ~ dnorm(mean = mu_rec[i],
                      sd = sigma_rec[ variable_rec[i] ])
  }
  
})

code_gaussian <- nimbleCode({
  
  ###############
  #  theta_inf  # 
  ###############
  
  mean_theta_inf ~ dnorm(0, sd = 10)
  sigma_theta_inf ~ T(dt(mu = 0, tau = 1, df = 3), 0, Inf)
  
  for(j in 1:n_var){
    beta_inf_raw[j] ~ dnorm(0, sd = 1)
    theta_inf[j] <- mean_theta_inf + beta_inf_raw[j] * sigma_theta_inf
    sigma_raw_old[j] ~ T(dt(mu = 0, tau = 1, df = 3), 0, Inf)
    sigma_old[j] <- s_old[j] * sigma_raw_old[j] 
    tau_old[j] <- 1/pow(sigma_old[j], 2)
  }
  
  ###############
  #  theta_0  # 
  ###############
  
  mean_theta_0 ~ dnorm(0, sd = 10)
  sigma_theta_0 ~ T(dt(mu = 0, tau = 1, df = 3), 0, Inf)
  
  for(j in 1:n_var){
    beta_0_raw[j] ~ dnorm(0, sd = 1)
    theta_0[j] <- mean_theta_0 + beta_0_raw[j] * sigma_theta_0
    sigma_raw_rec[j] ~ T(dt(mu = 0, tau = 1, df = 3), 0, Inf)
    sigma_rec[j] <- s_rec[j] * sigma_raw_rec[j] 
    tau_rec[j] <- 1/pow(sigma_rec[j], 2)
  }
  
  
  ############################
  # alpha_con, beta_con      #
  ############################
  
  mean_alpha_con  ~ dnorm(0, sd = 1)
  sigma_alpha_con ~ T(dt(mu = 0, tau = 1, df = 3), 0, Inf)
  
  mean_beta_con  ~ dnorm(0, sd = 1)
  sigma_beta_con ~ T(dt(mu = 0, tau = 1, df = 3), 0, Inf)
  
  
  for (j in 1:n_var) {
    alpha_con_raw[j] ~ dnorm(0, sd = 1)
    alpha_con[j] <- mean_alpha_con + alpha_con_raw[j] * sigma_alpha_con
    
    beta_con_raw[j] ~ dnorm(0, sd = 1)
    beta_con[j] <- mean_beta_con + beta_con_raw[j] * sigma_beta_con
  }
  
  ##########################
  # Likelihoods            #
  ##########################
  
  # old-growth
  for (i in 1:n_old) {
    Y_old[i] ~ dnorm(mean = theta_inf[ variable_old[i] ],
                      sd = sigma_old[ variable_old[i] ])
  }
  
  for (i in 1:n_rec) {
    
    mu_rec[i] <-
      theta_0[ variable_rec[i] ] +
      (theta_inf[ variable_rec[i] ] - theta_0[ variable_rec[i] ]) *
      (1 - exp(-lambda[i] * tx[i]))

    lambda[i] <- exp(
      alpha_con[ variable_rec[i] ] +
        beta_con[  variable_rec[i] ] * connectivity[i]
    )
    
    Y_rec[i] ~ dnorm(mean = mu_rec[i],
                      sdlog = sigma_rec[ variable_rec[i] ])
  }
  
})

# function to set initial values

initsFun_gaussian_strat <- function(constList, dataList) {
  n_var <- constList$n_var
  n_groups <- constList$n_groups # <--- IMPORTANTE: Pegar o número de grupos (3)
  
  list(
    # --- HIPERPARÂMETROS (Agora são vetores de tamanho n_groups) ---
    # Antes era apenas 0 ou runif(1), agora geramos um valor para cada grupo
    
    mean_theta_inf  = rnorm(n_groups, 0, 0.5), 
    mean_theta_0    = rnorm(n_groups, 0, 0.5),
    mean_alpha_con  = rnorm(n_groups, 0, 0.5),
    mean_beta_con   = rnorm(n_groups, 0, 0.5),
    
    sigma_theta_inf = rexp(n_groups, 1),
    sigma_theta_0   = rexp(n_groups, 1),
    sigma_alpha_con = rexp(n_groups, 1),
    sigma_beta_con  = rexp(n_groups, 1),
    
    # --- PARÂMETROS POR VARIÁVEL (Continuam iguais, tamanho n_var) ---
    beta_inf_raw  = rnorm(n_var, 0, 1),
    beta_0_raw    = rnorm(n_var, 0, 1),
    alpha_con_raw = rnorm(n_var, 0, 1),
    beta_con_raw  = rnorm(n_var, 0, 1),
    
    # --- ERROS DE OBSERVAÇÃO (Continuam iguais, tamanho n_var) ---
    sigma_raw_old = rexp(n_var, 1),
    sigma_raw_rec = rexp(n_var, 1)
  )
}


initsFun_gaussian <- function(constList, dataList) {
  n_var <- constList$n_var
  
  list(
    # hierarchical means
    mean_theta_inf  = 0,
    mean_theta_0    = 0,
    mean_alpha_con  = 0,
    mean_beta_con   = 0,
    
    # hierarchical sigmas (must be > 0)
    sigma_theta_inf = runif(1, 0.5, 1.5),
    sigma_theta_0   = runif(1, 0.5, 1.5),
    sigma_alpha_con = runif(1, 0.5, 1.5),
    sigma_beta_con  = runif(1, 0.5, 1.5),
    
    # group-level raw deviations (normal(0,1))
    beta_inf_raw  = rnorm(n_var, 0, 1),
    beta_0_raw    = rnorm(n_var, 0, 1),
    alpha_con_raw = rnorm(n_var, 0, 1),
    beta_con_raw  = rnorm(n_var, 0, 1),
    
    # standard deviation raw multipliers (half-t → must be positive)
    sigma_raw_old = runif(n_var, 0.5, 1.5),
    sigma_raw_rec = runif(n_var, 0.5, 1.5)
  )
}

monitorList <- c(
  "theta_0", "theta_inf",
  "alpha_con", "beta_con",
  "mean_theta_inf", "sigma_theta_inf",
  "mean_theta_0",   "sigma_theta_0",
  "mean_alpha_con", "sigma_alpha_con",
  "mean_beta_con",  "sigma_beta_con",
  "sigma_old", "sigma_rec"
)

# bundle the data for nimble
nimbleList <- list(
  code = code_gaussian_stratified,
  initsFun = initsFun_gaussian_strat,
  constList = constList,
  dataList = dataList,
  monitorList = monitorList
)

# set up model
modelR <- nimbleModel(
  code = code_gaussian_stratified,
  constants = constList,
  data = dataList,
  inits = initsFun_gaussian_strat(constList, dataList),
  calculate = FALSE
)
#modelR$initializeInfo()

# configure MCMC samplers & build MCMC
mcmcConf <-
  configureMCMC(modelR,
                monitors = monitorList,
                print = TRUE,
                nodes = NULL)
mcmcConf$replaceSampler(
  target = modelR$getNodeNames(
    stochOnly = TRUE,
    includeData = FALSE,
    includePredictive = FALSE
  ),
  type = "AF_slice"
)
mcmcConf$getUnsampledNodes()
mcmcR <- buildMCMC(mcmcConf)

# compile model, functions, set initial values and compile MCMC sampler
modelC <- compileNimble(modelR)
modelC$setInits(initsFun_gaussian_strat(constList, dataList))
mcmcC <- compileNimble(mcmcR, project = modelR)

# run the MCMC algorithm
samples_gaussian <- runMCMC(
  mcmc = mcmcC,
  niter = 2e4,
  nburnin = 1e2,
  thin = 3e1,
  nchains = 5,
  samplesAsCodaMCMC = TRUE
)

summary(samples_gaussian)

samples_gaussian_all <- as.mcmc(do.call(rbind, samples_gaussian))
colnames(samples_gaussian_all)

# diagnostics


gelman.diag(samples_gaussian[, grep("theta_inf", varnames(samples_gaussian)), drop=FALSE])
effectiveSize(samples_gaussian_all[, grepl("theta_inf", colnames(samples_gaussian_all)), drop = FALSE])
autocorr.diag(samples_gaussian_all[, grepl("theta_inf", colnames(samples_gaussian_all)), drop = FALSE])

effectiveSize(samples_gaussian_all[, grepl("theta_0", colnames(samples_gaussian_all)), drop = FALSE])
autocorr.diag(samples_gaussian_all[, grepl("theta_0", colnames(samples_gaussian_all)), drop = FALSE])

effectiveSize(samples_gaussian_all[, grepl("alpha_con", colnames(samples_gaussian_all)), drop = FALSE])
autocorr.diag(samples_gaussian_all[, grepl("alpha_con", colnames(samples_gaussian_all)), drop = FALSE])

effectiveSize(samples_gaussian_all[, grepl("beta_con", colnames(samples_gaussian_all)), drop = FALSE])
autocorr.diag(samples_gaussian_all[, grepl("beta_con", colnames(samples_gaussian_all)), drop = FALSE])

MCMCtrace(samples_gaussian, params = "theta_inf", ISB = F, exact = F, pdf = F)
MCMCtrace(samples_gaussian, params = "theta_0", ISB = F, exact = F, pdf = F)

MCMCtrace(samples_gaussian, params = "alpha_con", ISB = F, exact = F, pdf = F)
MCMCtrace(samples_gaussian, params = "beta_con", ISB = F, exact = F, pdf = F)

# plot
## probability distributions

# extract variables for plotting

sum_stats_g <- summary(samples_gaussian)$statistics[, "Mean"]

theta_0_est_g   <- sum_stats_g[grep("^theta_0\\[", names(sum_stats_g))]
theta_inf_est_g <- sum_stats_g[grep("^theta_inf\\[", names(sum_stats_g))]
alpha_con_est_g <- sum_stats_g[grep("^alpha_con\\[", names(sum_stats_g))]
beta_con_est_g  <- sum_stats_g[grep("^beta_con\\[", names(sum_stats_g))]

# raw data:

var_names_g <- colnames(dataSub)[-c(1:3)]

plot_data_g <- rec_df %>%
  mutate(VarName = factor(variable_rec, 
                          levels = 1:constList$n_var, 
                          labels = var_names))

rec_conn <- dataSub$ConIndex[dataSub$type == "rec"]
conn_values <- c(
  Low    = quantile(rec_conn, 0.25, na.rm = TRUE),
  Medium = quantile(rec_conn, 0.50, na.rm = TRUE),
  High   = quantile(rec_conn, 0.75, na.rm = TRUE)
)

#### Recovery time estimation ####

w <- weights_df %>% pull(prctg) # check if order is the same as in the model_df data frame

# extracting posterior samples for metrics

t_multi_conn_pol <- lapply(conn_values, function(cn) {
  recovery_tmulti_nimble_w(samples, groups = 1:3, conn = cn, weights = w[1:3])
})
sum(is.na(t_multi_conn_pol))

t_multi_conn_sd <- lapply(conn_values, function(cn) {
  recovery_tmulti_nimble_w(samples, groups = 4:6, conn = cn, weights = w[4:6])
})
sum(is.na(t_multi_conn_sd))

# per group:

t_per_group_conn_pol_g <- lapply(conn_values, function(cn) {
  recovery_tmulti_per_group_nimble(samples_gaussian, groups = 1:6, conn = cn)
})
sum(is.na(t_per_group_conn_pol_g[[1]]))
sum(is.na(t_per_group_conn_pol_g[[2]]))
sum(is.na(t_per_group_conn_pol_g[[3]]))
sum(is.na(t_per_group_conn_pol_g[[4]]))
sum(is.na(t_per_group_conn_pol_g[[5]]))
sum(is.na(t_per_group_conn_pol_g[[6]]))


t_per_group_conn_sd_g <- lapply(conn_values, function(cn) {
  recovery_tmulti_per_group_nimble(samples_gaussian, groups = 7:12, conn = cn)
})
sum(is.na(t_per_group_conn_sd[[1]]))
sum(is.na(t_per_group_conn_sd[[2]]))
sum(is.na(t_per_group_conn_sd[[3]]))

### credible intervals:

qt90_multi_conn_pol_g <- lapply(t_multi_conn_pol_g, function(x) {
  quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
})
qt90_multi_conn_pol

qt90_per_group_conn_pol_g <- lapply(t_per_group_conn_pol_g, function(tmat) {
  apply(tmat, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
        na.rm = TRUE)
})
qt90_per_group_conn_pol_g

qt90_multi_conn_sd_g <- lapply(t_multi_conn_sd_g, function(x) {
  quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
})
qt90_multi_conn_sd_g

qt90_per_group_conn_sd_g <- lapply(t_per_group_conn_sd_g, function(tmat) {
  apply(tmat, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
        na.rm = TRUE)
})
qt90_per_group_conn_sd_g

sigmaSq_rec <- do.call(rbind, as.mcmc.list(samp1$sigmaSq_rec))
sigmaSq_old <- do.call(rbind, as.mcmc.list(samp1$sigmaSq_old))


colnames(t90) <- as.vector(outer(names(metricList_conn), names(metricList_conn[[1]]), paste, sep = "_"))
qtheta_0 <- apply(theta_0, 2, FUN = quantile, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))
qtheta_inf <- apply(theta_inf, 2, FUN = quantile, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))
qalpha <- apply(alpha, 2, FUN = quantile, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))
qbeta <- apply(beta, 2, FUN = quantile, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))
qsigmarec <- apply(sigmaSq_rec, 2, FUN = quantile, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))
qsigmaold <- apply(sigmaSq_old, 2, FUN = quantile, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))

pvalue_beta <- c(mean(sign(beta[,1]) == sign(median(beta[,1]))),
                 mean(sign(beta[,2]) == sign(median(beta[,2]))),
                 mean(sign(beta[,3]) == sign(median(beta[,3]))))

# plants vs. animals
perc_pxa <- c(mean(metricList_conn$AlphaAnimals$Low >= metricList_conn$AlphaPlants$Low),
              mean(metricList_conn$AlphaAnimals$Medium >= metricList_conn$AlphaPlants$Medium),
              mean(metricList_conn$AlphaAnimals$High >= metricList_conn$AlphaPlants$High)) 
# plants vs. interactions
perc_pxi <- c(mean(metricList_conn$AlphaInt$Low >= metricList_conn$AlphaPlants$Low),
              mean(metricList_conn$AlphaInt$Medium >= metricList_conn$AlphaPlants$Medium),
              mean(metricList_conn$AlphaInt$High >= metricList_conn$AlphaPlants$High))
# animals vs. interactions
perc_axi <- c(mean(metricList_conn$AlphaInt$Low <= metricList_conn$AlphaAnimals$Low),
              mean(metricList_conn$AlphaInt$Medium <= metricList_conn$AlphaAnimals$Medium),
              mean(metricList_conn$AlphaInt$High <= metricList_conn$AlphaAnimals$High))

#### Save data ####

scores <- list(scores_int_hbt, scores_plants_hbt, scores_animals_hbt)
originalities <- list(alpha_int, alpha_plants, alpha_animals)
model_samples <- list(alpha, beta, theta_0, theta_inf, t90, qtheta_inf)
diff_tests <- list(qt90_conn, perc_pxa, perc_pxi, perc_axi, pvalue_beta)

# saveRDS(scores, "scores.RData")
# saveRDS(originalities, "originalities.RData")
# saveRDS(data1, "analysis_data.RData")
# saveRDS(model_samples, "model_samples.RData")
# saveRDS(diff_tests, "groups_diffs.RData")# saveRDS(scores, "scores.RData")
# saveRDS(originalities, "originalities.RData")
# saveRDS(data1, "analysis_data.RData")
# saveRDS(model_samples, "model_samples.RData")
# saveRDS(diff_tests, "groups_diffs.RData")


