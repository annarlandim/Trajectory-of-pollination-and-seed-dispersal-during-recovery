#### Data ####

data <- read.csv("data/processed/model_df.3.csv")

weights_df <- read.csv("data/processed/weights_df.csv")

data$type <- factor(ifelse(1:nrow(data) %in% grep("OG", data$Plot_ID), "old", "rec"),
                     levels = c("old", "rec"))

dataSub <- subset(data, select = c(type, RegTime, ConIndex,  
                                   FDBees, FDMoths,  FDBat_pol, 
                                   FDBats, FDBirds, FDNf#, 
                                   # FDSeeds, AbSeeds, RichSeeds, ShnSeeds,
                                   # FDSdlng, AbSdlng, RichSdlng, ShnSdlng
                                   )) 
table(dataSub$type)
str(dataSub)

# plot the data

# par(mfrow = c(1, 3))
# plot(FDBees ~ RegTime, dataSub)
# plot(FDMoths ~ RegTime, dataSub)
# plot(FDBat_pol ~ RegTime, dataSub)
# dev.off()
# 
# par(mfrow = c(1, 3))
# plot(FDBats ~ RegTime, dataSub)
# plot(FDBirds ~ RegTime, dataSub)
# plot(FDNf ~ RegTime, dataSub)
# dev.off()

# par(mfrow = c(1, 4))
# plot(AbSeeds ~ RegTime, dataSub)
# plot(RichSeeds ~ RegTime, dataSub)
# plot(ShnSeeds ~ RegTime, dataSub)
# plot(FDSeeds ~ RegTime, dataSub)
# dev.off()
# 
# par(mfrow = c(1, 4))
# plot(AbSdlng ~ RegTime, dataSub)
# plot(RichSdlng ~ RegTime, dataSub)
# plot(ShnSdlng ~ RegTime, dataSub)
# plot(FDSdlng ~ RegTime, dataSub)
# dev.off()

# define data for jags model

dataSub <- dataSub %>%
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

long <- dataSub %>%
  pivot_longer(cols = -c(1:3), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = colnames(dataSub)[-c(1:3)]), variable = as.integer(variable)) %>%
  filter(!is.na(value))

sd_old_vec <- long %>%
  filter(type == "old") %>%
  group_by(variable) %>%
  summarise(sdlog = sd(log(value), na.rm = TRUE), .groups = "drop") %>%
  pull(sdlog)

# to avoid very tiny sds:
which(sd_old_vec < 0.2)
sd_old_vec <- pmax(sd_old_vec, 0.2)

sd_rec_vec <- long %>%
  filter(type == "rec") %>%
  group_by(variable) %>%
  summarise(sdlog = sd(log(value), na.rm = TRUE), .groups = "drop") %>%
  pull(sdlog)

# to avoid very tiny sds:
which(sd_rec_vec < 0.2)
sd_rec_vec <- pmax(sd_rec_vec, 0.2)

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
  ) %>%
  mutate(connectivity = as.numeric(scale(connectivity)))
  

colnames(dataSub)[-c(1:3)]
group_index_vec <- c(1, 1, 1, 2, 2, 2)#, 3, 3, 3, 3, 4, 4, 4, 4) # 3 SD, 2Pol, 4Seeds, 4 Seedlings

constList <- list(
  variable_old = old_df$variable_old,
  variable_rec = rec_df$variable_rec,
  n_var = ncol(dataSub) - 3,
  n_groups = length(unique(group_index_vec)),  
  group_index = group_index_vec,
  n_old = nrow(old_df),
  n_rec = nrow(rec_df)
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

code_stratified <- nimbleCode({
  
  ################
  #  Hyperpriors # 
  ################
  
  for(k in 1:n_groups){
    
    # --- theta_inf ---
    mean_theta_inf[k]  ~ dnorm(0, sd = 2)
    sigma_theta_inf[k] ~  T(dt(0, pow(0.5, -2), 4), 0, Inf) 
    
    # --- theta_0 ---
    mean_theta_0[k]    ~ dnorm(0, sd = 2)
    sigma_theta_0[k]   ~  T(dt(0, pow(0.5, -2), 4), 0, Inf)
    
    # --- Connectivity (alpha and beta) ---
    mean_alpha_con[k]  ~ dnorm(-2.2, sd = 0.9) 
    sigma_alpha_con[k] ~ T(dt(0, pow(0.6, -2), 4), 0, Inf) 
    
    mean_beta_con[k]   ~ dnorm(0, sd = 0.35)
    sigma_beta_con[k]  ~ T(dt(0, pow(0.30, -2), 4), 0, Inf)
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
  n_groups <- constList$n_groups 
  
  list(
    
    mean_theta_inf  = rnorm(n_groups, 0, 0.1), 
    mean_theta_0    = rnorm(n_groups, 0, 0.1),
    mean_alpha_con  = rnorm(n_groups, -2.2, 0.1),
    # mean_alpha_con  = rnorm(n_groups, 0, 0.1),
    mean_beta_con   = rnorm(n_groups, 0, 0.1),
    
    # sigma_theta_inf = rexp(n_groups, 2),
    # sigma_theta_0   = rexp(n_groups, 2),
    # sigma_alpha_con = rexp(n_groups, 2),
    # sigma_beta_con  = rexp(n_groups, 2),
    
    sigma_theta_inf = rexp(n_groups, 1),
    sigma_theta_0   = rexp(n_groups, 1),
    sigma_alpha_con = rexp(n_groups, 1),
    sigma_beta_con  = rexp(n_groups, 1),
    
    beta_inf_raw  = rnorm(n_var, 0, 1),
    # beta_0_raw    = rnorm(n_var, 0, 0.5),
    beta_0_raw    = rnorm(n_var, 0, 1),
    alpha_con_raw = rnorm(n_var, 0, 1),
    # beta_con_raw  = rnorm(n_var, 0, 0.5),
    beta_con_raw  = rnorm(n_var, 0, 1),

    sigma_raw_old = rexp(n_var, 1),
    sigma_raw_rec = rexp(n_var, 1)
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

# Run parallel chains
# 1. Create a function that runs a single chain
run_one_chain <- function(seed, nimbleList) {
  library(nimble)
  # set.seed(seed)
  
  # Re-build the model inside the worker
  modelR <- nimbleModel(
    code = nimbleList$code,
    constants = nimbleList$constList,
    data = nimbleList$dataList,
    inits = nimbleList$initsFun(nimbleList$constList, nimbleList$dataList),
    calculate = FALSE
  )
  
  # Configure MCMC
  mcmcConf <- configureMCMC(modelR,
                            monitors = nimbleList$monitorList,
                            print = TRUE,
                            nodes = NULL)

  # use slice APENAS onde faz sentido
  # mcmcConf$replaceSampler(
  #   target = c("theta_0", "theta_inf", "alpha_con", "beta_con"),
  #   type = "AF_slice"
  # )
  # 
  # Apply your specific sampler choice
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
  
  modelC <- compileNimble(modelR)
  modelC$setInits(nimbleList$initsFun(nimbleList$constList, nimbleList$dataList))
  mcmcC <- compileNimble(mcmcR, project = modelR)
  
  # Run MCMC
  samples <- runMCMC(
    mcmc = mcmcC,
    niter = 1e5,   # Adjusted for efficiency
    nburnin = 2e4,
    thin = 5e1,        # Higher thinning helps reduce memory bottleneck
    samplesAsCodaMCMC = TRUE
  )
  
  return(samples)
}

this_cluster <- makeCluster(5)

clusterSetRNGStream(this_cluster, 1212)

# 3. Export your data and function to the cluster
clusterExport(this_cluster, c("nimbleList", "run_one_chain"))

# 4. Run in parallel
samples_raw <- parLapply(cl = this_cluster, X = 1:5, fun = run_one_chain, nimbleList = nimbleList)

# 5. Stop the cluster and convert to mcmc.list
stopCluster(this_cluster)
samples <- as.mcmc.list(samples_raw)

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

#### Recovery time estimation ####

w <- weights_df %>% pull(prctg) # check if order is the same as in the model_df data frame

rec_conn <- dataSub$ConIndex[dataSub$type == "rec"]
conn_values <- c(
  Low    = quantile(rec_conn, 0.25, na.rm = TRUE),
  Medium = quantile(rec_conn, 0.50, na.rm = TRUE),
  High   = quantile(rec_conn, 0.75, na.rm = TRUE)
)

### extracting posterior samples for metrics

t_multi_conn_pol <- lapply(conn_values, function(cn) {
  recovery_tmulti_nimble_w_con(samples, groups = 1:3, conn = cn, weights = w[1:3])
})
sum(is.na(t_multi_conn_pol))

t_multi_conn_sd <- lapply(conn_values, function(cn) {
  recovery_tmulti_nimble_w_con(samples, groups = 4:6, conn = cn, weights = w[4:6])
})
sum(is.na(t_multi_conn_sd))

# per group:

t_per_group_conn_pol <- lapply(conn_values, function(cn) {
  recovery_tmulti_per_group_nimble_con(samples, groups = 1:3, conn = cn)
})
sum(is.na(t_per_group_conn_pol[[1]]))
sum(is.na(t_per_group_conn_pol[[2]]))
sum(is.na(t_per_group_conn_pol[[3]]))


t_per_group_conn_sd <- lapply(conn_values, function(cn) {
  recovery_tmulti_per_group_nimble_con(samples, groups = 4:6, conn = cn)
})
sum(is.na(t_per_group_conn_sd[[1]]))
sum(is.na(t_per_group_conn_sd[[2]]))
sum(is.na(t_per_group_conn_sd[[3]]))

# seeds:

t_seeds_conn <- lapply(conn_values, function(cn) {
  recovery_tmulti_per_group_nimble_con(samples, groups = 7:10, conn = cn)
})
sum(is.na(t_seeds_conn[[1]]))
sum(is.na(t_seeds_conn[[2]]))
sum(is.na(t_seeds_conn[[3]]))

# seedlings:

t_seedlings_conn <- lapply(conn_values, function(cn) {
  recovery_tmulti_per_group_nimble_con(samples, groups = 11:14, conn = cn)
})
sum(is.na(t_seedlings_conn[[1]]))
sum(is.na(t_seedlings_conn[[2]]))
sum(is.na(t_seedlings_conn[[3]]))

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

qt90_seeds_conn <- lapply(t_seeds_conn, function(tmat) {
  apply(tmat, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
        na.rm = TRUE)
})
qt90_seeds_conn

qt90_seedlings_conn <- lapply(t_seedlings_conn, function(tmat) {
  apply(tmat, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
        na.rm = TRUE)
})
qt90_seedlings_conn

#### Save data ####

rec_results <- list(
  pollination = list(
    multi = qt90_multi_conn_pol,
    per_group = qt90_per_group_conn_pol
  ),
  dispersal = list(
    multi = qt90_multi_conn_sd,
    per_group = qt90_per_group_conn_sd
  )#,
  # plants = list(
  #   seeds = qt90_seeds_conn,
  #   seedlings = qt90_seedlings_conn
  # )
)
# saveRDS(rec_results, "output/t90.4_con.rds")
# saveRDS(samples, "output/model_posteriors.4_con.rds")
# metadata <- list(var_names = colnames(dataSub)[-c(1,3)],
#                  group_index = group_index_vec)
# saveRDS(metadata, "output/model_metadata_con.rds")


