library(dplyr)
library(tidyr)

library(rjags)
library(coda)
library(lattice)
load.module("glm")

#### Recovery time ####

data <- read.csv("data/processed/model_df.csv")

data <- data %>% mutate(ConIndex = scale(ConIndex))

data$type <- factor(ifelse(1:nrow(data) %in% grep("OG", data$Plot_ID), "old", "rec"),
                     levels = c("old", "rec"))

dataSub <- subset(data, select = c(type, RegTime, ConIndex, VerticalVH, MaxTH, AGB, FDBats, FC1Bats, FC2Bats, FDBirds, FC1Birds, FC2Birds, FDNf, FC1Nf, FC2Nf))
table(dataSub$type)
str(dataSub)

eps <- 1e-6
dataSub$VerticalVH <- pmax(dataSub$VerticalVH, eps)
dataSub$MaxTH      <- pmax(dataSub$MaxTH, eps)
dataSub$AGB        <- pmax(dataSub$AGB, eps)
dataSub$FDBats     <- pmax(dataSub$FDBats, eps)
dataSub$FDBirds    <- pmax(dataSub$FDBirds, eps)
dataSub$FDNf       <- pmax(dataSub$FDNf, eps)

# plot the data

par(mfrow = c(1, 3))
plot(VerticalVH ~ RegTime, dataSub)
plot(MaxTH ~ RegTime, dataSub)
plot(AGB ~ RegTime, dataSub)
dev.off()

par(mfrow = c(1, 3))
plot(FDBats ~ RegTime, dataSub)
plot(FDBirds ~ RegTime, dataSub)
plot(FDNf ~ RegTime, dataSub)
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
  select(-FC1Bats, -FC2Bats, -FC1Birds, -FC2Birds, -FC1Nf, -FC2Nf)

long <- dataSub %>%
  pivot_longer(cols = -c(1:3), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = colnames(dataSub)[-c(1:3)]), variable = as.integer(variable)) %>%
  filter(!is.na(value))

sd_old_vec <- long %>%
  filter(type == "old") %>%
  group_by(variable) %>%
  summarise(sdlog = sd(log(value), na.rm = TRUE), .groups = "drop") %>%
  pull(sdlog)

sd_rec_vec <- long %>%
  filter(type == "rec") %>%
  group_by(variable) %>%
  summarise(sdlog = sd(log(value), na.rm = TRUE), .groups = "drop") %>%
  pull(sdlog)

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

jagsData <- list(
  # old
  Y_old = old_df$Y_old,
  n_old = nrow(old_df),
  variable_old = old_df$variable_old,
  s_old = sd_old_vec,
  # rec
  Y_rec = rec_df$Y_rec,
  n_rec = nrow(rec_df),
  variable_rec = rec_df$variable_rec,
  s_rec = sd_rec_vec,
  tx = rec_df$tx,
  connectivity = as.numeric(rec_df$connectivity),
  # dimensões
  # had to change this manually because currently not working with functional composition:
  # n_var = ncol(dataSub) - 3 
  n_var = ncol(dataSub) - 3 
)


# jagsData <- 
#   with(dataSub, {
#     Y <- cbind(VerticalVH, MaxTH, AGB, FDBats, FDBirds, FDNf)
#     Y_old <- Y[type == "old", ]
#     Y_rec <- Y[type == "rec", ]
#     out <- list(
#       Y_old = Y_old,
#       n_old = nrow(Y_old),
#       s_old = apply(Y_old, 2, FUN = function(x) sd(log(x), na.rm = TRUE)),
#       Y_rec = Y_rec,
#       n_rec = nrow(Y_rec),
#       s_rec = apply(Y_rec, 2, FUN = function(x) sd(log(x), na.rm = TRUE)),
#       tx = RegTime[type == "rec"],
#       connectivity = dataSub$ConIndex[dataSub$type == "rec"],
#       nY = ncol(Y)
#     )
#     return(out)
#   })

# setup JAGS code for the model

# first new version:
# jagsModel <- 
#   "model{
#     # loop on the three different response variables
#     for (j in 1:nY) {
#       # Loop on observations in old-growth forests to calculate the predicted values
#       for (i in 1:n_old) {
#         # Model Likelihood for old-growth forests
#         Y_old[i, j] ~ dlnorm(log(theta_inf[j]), tau_old[j])
#       }
#       # Loop on observations in secondary forests to calculate the predicted values
#       for (i in 1:n_rec) {							
#         mu_rec[i, j] <-
#           theta_0[j] + 
#           (theta_inf[j] - theta_0[j]) *
#           (1 - exp(-lambda[i, j] * tx[i]))
#         # Model lambda depending on connectivity:
#         lambda[i,j] <- exp(alpha[j] + beta[j] * connectivity[i])
#         # Model Likelihood for secondary forests
#         Y_rec[i, j] ~ dlnorm(log(mu_rec[i, j]), tau_rec[j])
#       }
#       # priors on variance components
#       tau_old[j] ~ dscaled.gamma(s_old[j], 2)
#       sigmaSq_old[j] <- pow(tau_old[j], -1)
#       # sigma_old[j] ~ dunif(0, 5) # gpt suggestion. using sample sd does not work due to many NAs. very weakly informative.
#       # tau_old[j] <- 1 / pow(sigma_old[j], 2)
#       tau_rec[j] ~ dscaled.gamma(s_rec[j], 2)
#       sigmaSq_rec[j] <- pow(tau_rec[j], -1)
#       # sigma_rec[j] ~ dunif(0, 5)
#       # tau_rec[j] <- 1 / pow(sigma_rec[j], 2)
#       # prior on asymptotic attribute value (the two step approach below is )
#       log_theta_inf_raw[j] ~ dnorm(0, 1) # lognormal with sd of 1
#       # here is where the sigma_theta_inf should be defined, but not in every loop.
#       log_theta_0_raw[j] ~ dnorm(0, 1) # lognormal with sd of 1
#       # here is where the sigma_theta_0 should be defined, but not in every loop.
#       theta_inf[j] <- exp(log_theta_inf_raw[j] * sigma_theta_inf)
#       # prior on initial attribute value
#       theta_0[j] <- exp(log_theta_0_raw[j] * sigma_theta_0)
#       # prior on intercept and slope for lambda
#       alpha[j] ~ dnorm(0,1)
#       beta[j] ~ dnorm(0,1)
#     }
#     sigma_theta_inf ~ dscaled.gamma(1, 2) # re-scale it with sigma
#     sigma_theta_0 ~ dscaled.gamma(1, 2) # re-scale it with sigma
#   }"
# cat(jagsModel, file = "jagsModel.txt")


# function to calculate recovery time based on Poorter et al. (2021)
# T90 was calculated by calculating for each moment in time the absolute
# attribute value using the site-specific model equations
# Recovery time is defined as the time needed to recover to 90% of OGF values.

# -> see definition of tx in the function below

# recoveryFun <- function(samples, which = 1, conn, maxt = 1000) {
#   theta_0 <- do.call(rbind, as.mcmc.list(samples$theta_0))[, paste("theta_0", "[", which, "]", sep = "")]
#   theta_inf <- do.call(rbind, as.mcmc.list(samples$theta_inf))[, paste("theta_inf", "[", which, "]", sep = "")]
#   alpha <- do.call(rbind, as.mcmc.list(samples$alpha))[, paste("alpha", "[", which, "]", sep = "")]
#   beta  <- do.call(rbind, as.mcmc.list(samples$beta))[, paste("beta", "[", which, "]", sep = "")]
#   
#   # Calculate lambda based on connectivity
#   lambda <- exp(alpha + beta * conn)
#   
#   tx <- seq(1, maxt, 0.1)
#   
#   # Initialize vector to store t90 for each sample
#   t90_values <- numeric(length(theta_0))
#   
#   failed_indices <- c()
#   
#   # Loop through the samples to compute t90 for each one
#   for (x in 1:length(theta_0)) {
#     theta_t <- theta_0[x] + (theta_inf[x] - theta_0[x]) * (1 - exp(-lambda[x] * tx))
#     
#     if (theta_0[x] <= theta_inf[x]) {
#       t90_values[x] <- min(tx[which(theta_t > (0.9 * theta_inf[x]))])
#     } else {
#       t90_values[x] <- min(tx[which(theta_t < (1.1 * theta_inf[x]))])
#     }
#   }
#   
#   # Return the result as an MCMC object with just the t90 column
#   return(t90_values)
# }

# --- generic multifunctionality recovery time
# groups: integer indices of columns in Y to combine (e.g., 4:6 for SD; 1:3 for structure)
# conn:   scalar connectivity at which to evaluate t_multi (e.g., Q1/median/Q3)
# maxt/by/threshold: time grid and criterion

recovery_tmulti <- function(samples, groups, conn, maxt = 1000) {
  theta_0 <- do.call(rbind, as.mcmc.list(samples$theta_0))[ , paste0("theta_0[",   groups, "]"), drop = FALSE]
  theta_inf <- do.call(rbind, as.mcmc.list(samples$theta_inf))[, paste0("theta_inf[", groups, "]"), drop = FALSE]
  aMat <- do.call(rbind, as.mcmc.list(samples$alpha_con))[   , paste0("alpha_con[",     groups, "]"), drop = FALSE]
  bMat  <- do.call(rbind, as.mcmc.list(samples$beta_con))[   , paste0("beta_con[",     groups, "]"), drop = FALSE]
  
  tgrid  <- seq(0, maxt, 0.1)
  lambda <- exp(aMat + bMat * conn)         # [draw x group]
  t_multi <- rep(NA_real_, nrow(theta_0))       # one value per posterior draw
  
  for (tt in tgrid) {
    At  <- theta_0 + (theta_inf - theta_0) * (1 - exp(-lambda * tt))  # [draw x group]
    Mt  <- rowMeans(abs(At - theta_inf) / theta_inf)              # mean relative deviation
    hit <- is.na(t_multi) & (Mt <= 0.10)
    if (any(hit)) t_multi[hit] <- tt
    if (all(!is.na(t_multi))) break
  }
  t_multi
}

# function to set initial values
# jorg
initFun <- function(data) {
  out <- with(data, {
    list(
      alpha = runif(nY, 0, 1),
      beta = runif(nY, 0, 1)#,
      # theta_inf = runif(nY, 0, 1),
      # theta_0 = runif(nY, 0, 1),
      # tau_old = runif(nY, 0, 1),
      # tau_rec = runif(nY, 0, 1)
    )
  })
  return(out)
}

# new:
initFun <- function(jd) {
  function() {
    
    # 1e-3 avoid 0s in log
    
    # To start with theta in runif(0,1) as before I am setting its parameters in the
    # model to mean 0, sigma 1 and raw = log(init)
    
    list(
      # theta_inf = exp(mean + raw*sigma) -> escolha raw = log(init) com mean=0, sigma=1
      mean_theta_inf  = 0,
      sigma_theta_inf = 1,
      beta_inf_raw    = log(runif(jd$n_var, 1e-3, 1)),
      
      # theta_0 idem
      mean_theta_0  = 0,
      sigma_theta_0 = 1,
      beta_0_raw    = log(runif(jd$n_var, 1e-3, 1)),
      
      # alpha_con, beta_con começam ~ Unif(0,1) como antes
      mean_alpha_con  = 0,
      sigma_alpha_con = 1,
      alpha_con_raw   = runif(jd$n_var, 0, 1),
      
      mean_beta_con   = 0,
      sigma_beta_con  = 1,
      beta_con_raw    = runif(jd$n_var, 0, 1),
      
      # precisões positivas (opcional iniciar)
      tau_old = 1/pmax(jd$s_old^2, 1e-3),
      tau_rec = 1/pmax(jd$s_rec^2, 1e-3)
    )
  }
}

# without hyper prior:

initFun <- function(jd) {
  function() {
    
    # 1e-3 avoid 0s in log
    
    # To start with theta in runif(0,1) as before I am setting its parameters in the
    # model to mean 0, sigma 1 and raw = log(init)
    
    list(
      theta_inf  = runif(jd$n_var, 1e-3, 1),
      theta_0  = runif(jd$n_var, 1e-3, 1),
      
      # precisões positivas (opcional iniciar)
      tau_old = 1/pmax(jd$s_old^2, 1e-3),
      tau_rec = 1/pmax(jd$s_rec^2, 1e-3),
      
      # alpha_con, beta_con começam ~ Unif(0,1) como antes
      alpha_con  = runif(jd$n_var, 0, 1),
      beta_con    = runif(jd$n_var, 0, 1)
    )
  }
}

# initialize the model
## i.e., get probability distributions
set.seed(1212) 
init1 <- jags.model("jagsModel.txt",
                    data = jagsData,
                    inits = initFun(jagsData),
                    n.chains = 5, n.adapt = 1e4)

# sample from posterior distribution

#update(init1, n.iter = 1e4)
monitor1 <- c("theta_0", "theta_inf", "alpha", "beta", "sigmaSq_old", "sigmaSq_rec", "sigma_theta_inf", "sigma_theta_0")
monitor2 <- c("theta_0", "theta_inf", "alpha_con", "beta_con", "sigmaSq_old", "sigmaSq_rec", 
              "mean_theta_0", "mean_theta_inf", "sigma_theta_0", "sigma_theta_inf", "mean_alpha_con", "mean_beta_con",
              "sigma_alpha_con", "sigma_beta_con")
monitor3 <- c("theta_0", "theta_inf", "alpha_con", "beta_con", "tau_old", "tau_rec", "sigmaSq_old", "sigmaSq_rec")

samp1 <- jags.samples(init1, variable.names = monitor2, n.iter = 2e5, thin = 2e2)
## thin reduces autocorr between consecutive samples in MCMC
## in that way, although we run 4e4 iterations, we only save 4e4 / 1e2 * 5 samples

# diagnostics

gelman.diag(as.mcmc.list(log(samp1$theta_0)))
effectiveSize(as.mcmc.list(log(samp1$theta_0)))
autocorr.diag(as.mcmc.list(log(samp1$theta_0)))

gelman.diag(as.mcmc.list(log(samp1$theta_inf)))
effectiveSize(as.mcmc.list(log(samp1$theta_inf)))
autocorr.diag(as.mcmc.list(log(samp1$theta_inf)))

gelman.diag(as.mcmc.list((samp1$alpha_con)))
effectiveSize(as.mcmc.list((samp1$alpha_con)))
autocorr.diag(as.mcmc.list((samp1$alpha_con)))

gelman.diag(as.mcmc.list((samp1$beta_con)))
effectiveSize(as.mcmc.list((samp1$beta_con)))
autocorr.diag(as.mcmc.list((samp1$beta_con)))

gelman.diag(as.mcmc.list(log(samp1$sigmaSq_old)))
effectiveSize(as.mcmc.list(log(samp1$sigmaSq_old)))
autocorr.diag(as.mcmc.list(log(samp1$sigmaSq_old)))

gelman.diag(as.mcmc.list(log(samp1$sigmaSq_rec)))
effectiveSize(as.mcmc.list(log(samp1$sigmaSq_rec)))
autocorr.diag(as.mcmc.list(log(samp1$sigmaSq_rec)))

gelman.diag(as.mcmc.list((samp1$mean_theta_0)))
effectiveSize(as.mcmc.list((samp1$mean_theta_0)))
autocorr.diag(as.mcmc.list((samp1$mean_theta_0)))

gelman.diag(as.mcmc.list((samp1$mean_theta_inf)))
effectiveSize(as.mcmc.list((samp1$mean_theta_inf)))
autocorr.diag(as.mcmc.list((samp1$mean_theta_inf)))
traceplot(as.mcmc.list(samp1$sigma_theta_inf))

gelman.diag(as.mcmc.list(log(samp1$sigma_theta_0)))
effectiveSize(as.mcmc.list(log(samp1$sigma_theta_0)))
autocorr.diag(as.mcmc.list(log(samp1$sigma_theta_0)))

gelman.diag(as.mcmc.list(log(samp1$sigma_theta_inf)))
effectiveSize(as.mcmc.list(log(samp1$sigma_theta_inf)))
autocorr.diag(as.mcmc.list(log(samp1$sigma_theta_inf)))

gelman.diag(as.mcmc.list((samp1$mean_alpha_con)))
effectiveSize(as.mcmc.list((samp1$mean_alpha_con)))
autocorr.diag(as.mcmc.list((samp1$mean_alpha_con)))

gelman.diag(as.mcmc.list((samp1$mean_beta_con)))
effectiveSize(as.mcmc.list((samp1$mean_beta_con)))
autocorr.diag(as.mcmc.list((samp1$mean_beta_con)))

gelman.diag(as.mcmc.list(log(samp1$sigma_alpha_con)))
effectiveSize(as.mcmc.list(log(samp1$sigma_alpha_con)))
autocorr.diag(as.mcmc.list(log(samp1$sigma_alpha_con)))

gelman.diag(as.mcmc.list(log(samp1$sigma_beta_con)))
effectiveSize(as.mcmc.list(log(samp1$sigma_beta_con)))
autocorr.diag(as.mcmc.list(log(samp1$sigma_beta_con)))

# plot
## probability distributions

densityplot(as.mcmc.list(samp1$alpha_con))
densityplot(as.mcmc.list(samp1$beta_con))
densityplot(as.mcmc.list(samp1$theta_0), scale = list(x = list(log = 10)))
densityplot(as.mcmc.list(samp1$theta_inf), scale = list(x = list(log = 10)))
densityplot(as.mcmc.list(samp1$sigmaSq_old), scale = list(x = list(log = 10)))
densityplot(as.mcmc.list(samp1$sigmaSq_rec), scale = list(x = list(log = 10)))
densityplot(as.mcmc.list(samp1$mean_theta_0))
densityplot(as.mcmc.list(samp1$mean_theta_inf))
densityplot(as.mcmc.list(samp1$sigma_theta_0), scale = list(x = list(log = 10)))
densityplot(as.mcmc.list(samp1$sigma_theta_inf), scale = list(x = list(log = 10)))
densityplot(as.mcmc.list(samp1$mean_alpha_con))
densityplot(as.mcmc.list(samp1$mean_beta_con))
densityplot(as.mcmc.list(samp1$sigma_alpha_con), scale = list(x = list(log = 10)))
densityplot(as.mcmc.list(samp1$sigma_beta_con), scale = list(x = list(log = 10)))


# extracting posterior samples for metrics

conn_values <- c(
  Low      = summary(jagsData$connectivity)[[2]],
  Medium   = summary(jagsData$connectivity)[[4]],
  High     = summary(jagsData$connectivity)[[5]]
)

metricList_conn <- lapply(1:3, function(metric) {
  lapply(conn_values, function(conn) recoveryFun(samp1, which = metric, conn = conn))
})
names(metricList_conn) <- colnames(jagsData$Y_old)

# extract variables for plotting

theta_inf <- do.call(rbind, as.mcmc.list(samp1$theta_inf))
theta_0 <- do.call(rbind, as.mcmc.list(samp1$theta_0))
alpha <- do.call(rbind, as.mcmc.list(samp1$alpha))
beta <- do.call(rbind, as.mcmc.list(samp1$beta))
t90 <- do.call(cbind, lapply(1:3, function(metric) {
  sapply(1:length(names(metricList_conn[[1]])), function(i) {
    as.numeric(metricList_conn[[metric]][[i]])
  })
}))
sigmaSq_rec <- do.call(rbind, as.mcmc.list(samp1$sigmaSq_rec))
sigmaSq_old <- do.call(rbind, as.mcmc.list(samp1$sigmaSq_old))

# qt90_conn <- lapply(metricList_conn, function(metric_data) {
#   lapply(metric_data, function(x) {
#     quantile(as.numeric(x), probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
#   })
# })

qtMulti_conn <- lapply(tmultiList_conn, function(proc_list) {
  lapply(proc_list, summ)
})


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


