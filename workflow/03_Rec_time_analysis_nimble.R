#### Data ####

data <- read.csv("data/processed/model_df.csv")

weights_df <- read.csv("data/processed/weights_df.csv")

data <- data %>% mutate(ConIndex = as.numeric(scale(ConIndex)))

data$type <- factor(ifelse(1:nrow(data) %in% grep("OG", data$Plot_ID), "old", "rec"),
                     levels = c("old", "rec"))

dataSub <- subset(data, select = c(type, RegTime, ConIndex,  
                                   FDBees, FC1Bees, FC2Bees, FDMoths, FC1Moths, FC2Moths, FDBat_pol, FC1Bat_pol, FC2Bat_pol,
                                   FDBats, FC1Bats, FC2Bats, FDBirds, FC1Birds, FC2Birds, FDNf, FC1Nf, FC2Nf,
                                   FDSdlng, FC1Sdlng, FC2Sdlng, AbSdlng, RichSdlng, ShnSdlng,
                                   VerticalVH, MaxTH, AGB)) 
table(dataSub$type)
str(dataSub)

# plot the data

par(mfrow = c(1, 3))
plot(FDBees ~ RegTime, dataSub)
plot(FDMoths ~ RegTime, dataSub)
plot(FDBat_pol ~ RegTime, dataSub)
dev.off()

par(mfrow = c(1, 3))
plot(FDBats ~ RegTime, dataSub)
plot(FDBirds ~ RegTime, dataSub)
plot(FDNf ~ RegTime, dataSub)
dev.off()

par(mfrow = c(1, 3))
plot(VerticalVH ~ RegTime, dataSub)
plot(MaxTH ~ RegTime, dataSub)
plot(AGB ~ RegTime, dataSub)
dev.off()

par(mfrow = c(1, 4))
plot(AbSdlng ~ RegTime, dataSub)
plot(RichSdlng ~ RegTime, dataSub)
plot(ShnSdlng ~ RegTime, dataSub)
plot(FDSdlng ~ RegTime, dataSub)
dev.off()

# define data for jags model

dataSub <- dataSub %>%
  select(-FC1Bees, -FC2Bees, -FC1Moths, -FC2Moths, -FC1Bat_pol, -FC2Bat_pol, -FC1Bats, -FC2Bats, -FC1Birds, -FC2Birds, -FC1Nf, -FC2Nf, -FC1Sdlng, -FC2Sdlng) %>%
  mutate(across(
    .cols = -c(1:3),     
    # scaling by dividing the max per column.
    # this is the best for logscale, where negatives are not allowed
    .fns = ~ .x / max(.x, na.rm = TRUE) 
  ))

str(dataSub)

zeros <- dataSub %>% 
  map_lgl(~ any(. == 0, na.rm = TRUE))
print(zeros)

# to avoid problems with log(0)
eps <- 1e-6
dataSub$AbSdlng <- pmax(dataSub$AbSdlng, eps)
dataSub$RichSdlng <- pmax(dataSub$RichSdlng, eps)
dataSub$ShnSdlng <- pmax(dataSub$ShnSdlng, eps)
dataSub$VerticalVH <- pmax(dataSub$VerticalVH, eps)
dataSub$MaxTH      <- pmax(dataSub$MaxTH, eps)
dataSub$AGB        <- pmax(dataSub$AGB, eps)

long <- dataSub %>%
  pivot_longer(cols = -c(1:3), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = colnames(dataSub)[-c(1:3)]), variable = as.integer(variable)) %>%
  filter(!is.na(value))

# mu_old_vec <- long %>%
#   filter(type == "old") %>%
#   group_by(variable) %>%
#   summarise(mulog = median(log(value), na.rm = TRUE), .groups = "drop") %>%
#   pull(mulog)
# 
sd_old_vec <- long %>%
  filter(type == "old") %>%
  group_by(variable) %>%
  summarise(sdlog = sd(log(value), na.rm = TRUE), .groups = "drop") %>%
  pull(sdlog)

# to avoid very tiny sds:
sd_old_vec <- pmax(sd_old_vec, 0.2)

# mu_0_vec <- long %>%
#   filter(RegTime == 0) %>%
#   group_by(variable) %>%
#   summarise(mu0 = median(log(value), na.rm = TRUE), .groups = "drop") %>%
#   pull(mu0)
# 
# sd_0_vec <- long %>%
#   filter(RegTime == 0) %>%
#   group_by(variable) %>%
#   summarise(s0 = sd(log(value), na.rm = TRUE), .groups = "drop") %>%
#   pull(s0)

sd_rec_vec <- long %>%
  filter(type == "rec") %>%
  group_by(variable) %>%
  summarise(sdlog = sd(log(value), na.rm = TRUE), .groups = "drop") %>%
  pull(sdlog)

sd_rec_vec <- pmax(sd_rec_vec, 0.2)

# Not enough data to estimate sd for non-flying mammals (n=1)
# sd_0_vec[!is.finite(sd_0_vec) ] <- sd_rec_vec[!is.finite(sd_0_vec) ]

# sd_0_vec <- pmax(sd_0_vec, 0.2)

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

group_index_vec <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4)

constList <- list(
  variable_old = old_df$variable_old,
  variable_rec = rec_df$variable_rec,
  n_var = ncol(dataSub) - 3,
  n_groups = 4,          # mudar quando add as plantulas
  group_index = group_index_vec,
  n_old = nrow(old_df),
  n_rec = nrow(rec_df)
)

dataList <- list(
  Y_old = old_df$Y_old,
  # mu_old = mu_old_vec,
  s_old = sd_old_vec,
  # mu_0 = mu_0_vec,
  # s_0 = sd_0_vec,
  Y_rec = rec_df$Y_rec,
  s_rec = sd_rec_vec,
  tx = rec_df$tx,
  connectivity = as.numeric(rec_df$connectivity)
)


#### Model ####

code_stratified <- nimbleCode({
  
  ################
  #  Hyperpriors # 
  ################
  
  for(k in 1:n_groups){
    
    # --- theta_inf ---
    mean_theta_inf[k]  ~ dnorm(0, sd = 2)
    sigma_theta_inf[k] ~ dexp(0.5)
    
    # --- theta_0 ---
    mean_theta_0[k]    ~ dnorm(0, sd = 2)
    sigma_theta_0[k]   ~ dexp(0.5)
    
    # --- Connectivity (alpha and beta) ---
    mean_alpha_con[k]  ~ dnorm(0, sd = 2)
    sigma_alpha_con[k] ~ T(dt(0, tau = 1, df = 3), 0, Inf)
    
    mean_beta_con[k]   ~ dnorm(0, sd = 2)
    sigma_beta_con[k]  ~  T(dt(0, tau = 1, df = 3), 0, Inf)
  }
  
  ################################
  #   Parameters per variable j  #
  ################################
  
  # --- theta_inf ---
  
  for(j in 1:n_var){
    beta_inf_raw[j] ~ dnorm(0, sd = 1)
    # asymptote per group, original scale (always > 0)
    theta_inf[j] <- exp(
      mean_theta_inf[ group_index[j] ] + 
        beta_inf_raw[j] * sigma_theta_inf[ group_index[j] ]
    )
    # priors on variance components:
    sigma_raw_old[j] ~ dexp(0.5)
    sigma_old[j] <- s_old[j] * sigma_raw_old[j] # scale by observed sd
    tau_old[j] <- 1/pow(sigma_old[j], 2)
  }
  
  # --- theta_0 ---
  
  for(j in 1:n_var){
    beta_0_raw[j] ~ dnorm(0, sd = 1)
    # asymptote per group, original scale (always > 0)
    theta_0[j] <- exp(
      mean_theta_0[ group_index[j] ] + 
        beta_0_raw[j] * sigma_theta_0[ group_index[j] ]
    )
    # priors on variance components:
    sigma_raw_rec[j] ~ dexp(0.5)
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
    Y_old[i] ~ dlnorm(meanlog = log(theta_inf[ variable_old[i] ]),
                      sdlog   = sigma_old[ variable_old[i] ])
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
    
    Y_rec[i] ~ dlnorm(meanlog = log(mu_rec[i]),
                      sdlog   = sigma_rec[ variable_rec[i] ])
  }
  
})

code <- nimbleCode({
  
  ###############
  #  theta_inf  # 
  ###############
  
  mean_theta_inf ~ dnorm(0, sd = 1)
  sigma_theta_inf ~ T(dt(mu = 0, tau = 1, df = 3), 0, Inf)
  
  for(j in 1:n_var){
    beta_inf_raw[j] ~ dnorm(0, sd = 1)
    # asymptote per group, original scale (always > 0)
    theta_inf[j] <- exp(mean_theta_inf + beta_inf_raw[j] * sigma_theta_inf)
    # priors on variance components:
    sigma_raw_old[j] ~ T(dt(mu = 0, tau = 1, df = 3), 0, Inf)
    sigma_old[j] <- s_old[j] * sigma_raw_old[j] # scale by observed sd
    tau_old[j] <- 1/pow(sigma_old[j], 2)
  }
  
  ###############
  #  theta_0  # 
  ###############
  
  mean_theta_0 ~ dnorm(0, sd = 1)
  sigma_theta_0 ~ T(dt(mu = 0, tau = 1, df = 3), 0, Inf)
  
  for(j in 1:n_var){
    beta_0_raw[j] ~ dnorm(0, sd = 1)
    # asymptote per group, original scale (always > 0)
    theta_0[j] <- exp(mean_theta_0 + beta_0_raw[j] * sigma_theta_0)
    # priors on variance components:
    sigma_raw_rec[j] ~ T(dt(mu = 0, tau = 1, df = 3), 0, Inf)
    sigma_rec[j] <- s_rec[j] * sigma_raw_rec[j] # scale by observed sd
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
    # parameterization reminder:
    # meanlog = log(theta_inf[group]), sdlog = sigma_old[group]
    Y_old[i] ~ dlnorm(meanlog = log(theta_inf[ variable_old[i] ]),
                      sdlog   = sigma_old[ variable_old[i] ])
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
    
    Y_rec[i] ~ dlnorm(meanlog = log(mu_rec[i]),
                      sdlog   = sigma_rec[ variable_rec[i] ])
  }
  
})

# function to set initial values

initsFun_strat <- function(constList, dataList) {
  n_var <- constList$n_var
  n_groups <- constList$n_groups # <--- IMPORTANTE: Pegar o número de grupos (3)
  
  list(
    # --- HIPERPARÂMETROS (Agora são vetores de tamanho n_groups) ---
    # Antes era apenas 0 ou runif(1), agora geramos um valor para cada grupo
    
    mean_theta_inf  = rnorm(n_groups, 0, 0.1), 
    mean_theta_0    = rnorm(n_groups, 0, 0.1),
    mean_alpha_con  = rnorm(n_groups, 0, 0.1),
    mean_beta_con   = rnorm(n_groups, 0, 0.1),
    
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

initsFun <- function(constList, dataList) {
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
  code = code_stratified,
  initsFun = initsFun_strat,
  constList = constList,
  dataList = dataList,
  monitorList = monitorList
)

# set up model
modelR <- nimbleModel(
  code = code_stratified,
  constants = constList,
  data = dataList,
  inits = initsFun_strat(constList, dataList),
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
modelC$setInits(initsFun_strat(constList, dataList))
mcmcC <- compileNimble(mcmcR, project = modelR)

# run the MCMC algorithm
samples <- runMCMC(
  mcmc = mcmcC,
  niter = 5e4,
  nburnin = 1e3,
  thin = 1e2,
  nchains = 5,
  samplesAsCodaMCMC = TRUE
)

summary(samples)

# use the median of the posterior as the intivalues. doesnt have to be a funciton, but I need to se how it is 
samples_all <- as.mcmc(do.call(rbind, samples))
colnames(samples_all)

# diagnostics

gelman_log <- function(samples_list, pattern) {

  names_vars <- varnames(samples_list)
  idx <- grep(paste0(pattern, "\\["), names_vars) 
  
  samples_log <- as.mcmc.list(lapply(samples_list, function(chain) {
    
    subset_chain <- chain[, idx, drop = FALSE]
    
    return(as.mcmc(log(subset_chain))) 
  }))
  return(gelman.diag(samples_log))
}

gelman_log(samples, "^theta_inf")
effectiveSize(log(samples_all[, grepl("^theta_inf\\[", colnames(samples_all)), drop = FALSE]))
autocorr.diag(log(samples_all[, grepl("^theta_inf\\[", colnames(samples_all)), drop = FALSE]))

gelman.diag(samples[, grep("^(mean_theta_inf|sigma_theta_inf)\\[", varnames(samples)), drop=FALSE])
effectiveSize(samples_all[, grep("^(mean_theta_inf|sigma_theta_inf)\\[", varnames(samples)), drop = FALSE])
autocorr.diag(samples_all[, grep("^(mean_theta_inf|sigma_theta_inf)\\[", varnames(samples)), drop = FALSE])

gelman_log(samples, "^theta_0")
effectiveSize(log(samples_all[, grepl("^theta_0\\[", colnames(samples_all)), drop = FALSE]))
autocorr.diag(log(samples_all[, grepl("^theta_0\\[", colnames(samples_all)), drop = FALSE]))

gelman.diag((samples[, grep("^(mean_theta_0|sigma_theta_0)\\[", varnames(samples)), drop=FALSE]))
effectiveSize(samples_all[, grep("^(mean_theta_0|sigma_theta_0)\\[", varnames(samples)), drop = FALSE])
autocorr.diag(samples_all[, grep("^(mean_theta_0|sigma_theta_0)\\[", varnames(samples)), drop = FALSE])

gelman.diag((samples[, grep("^alpha_con\\[", varnames(samples)), drop=FALSE]))
effectiveSize(samples_all[, grepl("alpha_con", colnames(samples_all)), drop = FALSE])
autocorr.diag(samples_all[, grepl("alpha_con", colnames(samples_all)), drop = FALSE])

gelman.diag((samples[, grep("^(mean_alpha_con|sigma_alpha_con)\\[", varnames(samples)), drop=FALSE]))
effectiveSize(samples_all[, grep("^(mean_alpha_con|sigma_alpha_con)\\[", varnames(samples)), drop = FALSE])
autocorr.diag(samples_all[, grep("^(mean_alpha_con|sigma_alpha_con)\\[", varnames(samples)), drop = FALSE])

gelman.diag((samples[, grep("^beta_con\\[", varnames(samples)), drop=FALSE]))
effectiveSize(samples_all[, grepl("beta_con", colnames(samples_all)), drop = FALSE])
autocorr.diag(samples_all[, grepl("beta_con", colnames(samples_all)), drop = FALSE])

gelman.diag((samples[, grep("^(mean_beta_con|sigma_beta_con)\\[", varnames(samples)), drop=FALSE]))
effectiveSize(samples_all[, grep("^(mean_beta_con|sigma_beta_con)\\[", varnames(samples)), drop = FALSE])
autocorr.diag(samples_all[, grep("^(mean_beta_con|sigma_beta_con)\\[", varnames(samples)), drop = FALSE])


MCMCtrace(samples, params = "theta_inf", ISB = F, exact = F, pdf = F)
MCMCtrace(samples, params = "theta_0", ISB = F, exact = F, pdf = F)

MCMCtrace(samples, params = "alpha_con", ISB = F, exact = F, pdf = F)
MCMCtrace(samples, params = "beta_con", ISB = F, exact = F, pdf = F)

##### plot

# extract variables for plotting

sum_stats <- summary(samples)$statistics[, "Mean"]

theta_0_est   <- sum_stats[grep("^theta_0\\[", names(sum_stats))]
theta_inf_est <- sum_stats[grep("^theta_inf\\[", names(sum_stats))]
alpha_con_est <- sum_stats[grep("^alpha_con\\[", names(sum_stats))]
beta_con_est  <- sum_stats[grep("^beta_con\\[", names(sum_stats))]

# raw data:

var_names <- colnames(dataSub)[-c(1:3)]

plot_data <- rec_df %>%
  mutate(VarName = factor(variable_rec, 
                          levels = 1:constList$n_var, 
                          labels = var_names))

rec_conn <- dataSub$ConIndex[dataSub$type == "rec"]
conn_values <- c(
  Low    = quantile(rec_conn, 0.25, na.rm = TRUE),
  Medium = quantile(rec_conn, 0.50, na.rm = TRUE),
  High   = quantile(rec_conn, 0.75, na.rm = TRUE)
)

# Create an empty list to store curve data
curve_list <- list()

for(j in 1:9) {
  
  # Parameters for this specific variable
  t0   <- unname(theta_0_est[j])
  tinf <- unname(theta_inf_est[j])
  a    <- unname(alpha_con_est[j])
  b    <- unname(beta_con_est[j])
  
  # Calculate Lambda 
  lambda  <- exp(a + b * conn_values[2]) 
  
  # Calculate Predicted Y
  y_pred <- t0 + (tinf - t0) * (1 - exp(-lambda * time_seq))
  
  curve_list[[j]] <- data.frame(
    VarName = var_names[j],
    tx = time_seq,
    Y_pred = y_pred,
    theta_inf = tinf 
  )
}

# Combine all curves into one dataframe
curve_df <- do.call(rbind, curve_list)
curve_df$VarName <- factor(curve_df$VarName, levels = var_names)

ggplot() +
  # 1. Plot the Raw Data Points
  geom_point(data = plot_data, aes(x = tx, y = Y_rec, color = connectivity), alpha = 0.6) +
  
  # 2. Plot the Fitted Curve (Mean Connectivity)
  geom_line(data = curve_df, aes(x = tx, y = Y_pred), size = 1, color = "black") +
  
  # 3. Plot the Asymptote (Optional but helpful)
  geom_hline(data = distinct(curve_df, VarName, theta_inf), 
             aes(yintercept = theta_inf), linetype = "dashed", color = "red") +
  
  # Aesthetics
  scale_color_gradient2(low = "red", mid = "gray", high = "blue", midpoint = 0,
                        name = "Connectivity\n(Scaled)") +
  facet_wrap(~VarName, scales = "free_y") +
  theme_bw() +
  labs(x = "Recovery Time (years)", 
       y = "Functional Diversity",
       title = "Model Fit: Predicted Recovery vs. Observed Data",
       subtitle = "Black line represents recovery at Selected Connectivity Level")

#### Recovery time estimation ####

w <- weights_df %>% pull(prctg) # check if order is the same as in the model_df data frame

### extracting posterior samples for metrics

t_multi_conn_pol <- lapply(conn_values, function(cn) {
  recovery_tmulti_nimble_w(samples, groups = 1:3, conn = cn, weights = w[1:3])
})
sum(is.na(t_multi_conn_pol))

t_multi_conn_sd <- lapply(conn_values, function(cn) {
  recovery_tmulti_nimble_w(samples, groups = 4:6, conn = cn, weights = w[4:6])
})
sum(is.na(t_multi_conn_sd))

# per group:

t_per_group_conn_pol <- lapply(conn_values, function(cn) {
  recovery_tmulti_per_group_nimble(samples, groups = 1:3, conn = cn)
})
sum(is.na(t_per_group_conn_pol[[1]]))
sum(is.na(t_per_group_conn_pol[[2]]))
sum(is.na(t_per_group_conn_pol[[3]]))


t_per_group_conn_sd <- lapply(conn_values, function(cn) {
  recovery_tmulti_per_group_nimble(samples, groups = 4:6, conn = cn)
})
sum(is.na(t_per_group_conn_sd[[1]]))
sum(is.na(t_per_group_conn_sd[[2]]))
sum(is.na(t_per_group_conn_sd[[3]]))

# seedlings:

t_seedlings_conn <- lapply(conn_values, function(cn) {
  recovery_tmulti_per_group_nimble(samples, groups = 7:10, conn = cn)
})
sum(is.na(t_seedlings_conn[[1]]))
sum(is.na(t_seedlings_conn[[2]]))
sum(is.na(t_seedlings_conn[[3]]))

# forest structure:

t_f_structure_conn <- lapply(conn_values, function(cn) {
  recovery_tmulti_per_group_nimble(samples, groups = 11:13, conn = cn)
})
sum(is.na(t_f_structure_conn[[1]]))
sum(is.na(t_f_structure_conn[[2]]))
sum(is.na(t_f_structure_conn[[3]]))

### credible intervals:

qt90_multi_conn_pol <- lapply(t_multi_conn_pol, function(x) {
  quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
})
qt90_multi_conn_pol

qt90_per_group_conn_pol <- lapply(t_per_group_conn_pol, function(tmat) {
  apply(tmat, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
    na.rm = TRUE)
})
qt90_per_group_conn_pol

qt90_multi_conn_sd <- lapply(t_multi_conn_sd, function(x) {
  quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
})
qt90_multi_conn_sd

qt90_per_group_conn_sd <- lapply(t_per_group_conn_sd, function(tmat) {
  apply(tmat, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
        na.rm = TRUE)
})
qt90_per_group_conn_sd

qt90_seedlings_conn <- lapply(t_seedlings_conn, function(tmat) {
  apply(tmat, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
        na.rm = TRUE)
})
qt90_seedlings_conn

qt90_f_structure_conn <- lapply(t_f_structure_conn, function(tmat) {
  apply(tmat, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
        na.rm = TRUE)
})
qt90_f_structure_conn


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


