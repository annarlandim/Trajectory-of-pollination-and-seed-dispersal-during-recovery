theta_0 <- do.call(rbind, as.mcmc.list(samp1$theta_0))[ , paste0("theta_0[",   4:6, "]"), drop = FALSE]
theta_inf <- do.call(rbind, as.mcmc.list(samp1$theta_inf))[, paste0("theta_inf[", 4:6, "]"), drop = FALSE]
aMat <- do.call(rbind, as.mcmc.list(samp1$alpha))[   , paste0("alpha[",     4:6, "]"), drop = FALSE]
bMat  <- do.call(rbind, as.mcmc.list(samp1$beta))[   , paste0("beta[",     4:6, "]"), drop = FALSE]

tgrid  <- seq(0, 1000, 0.1)
lambda <- exp(aMat + bMat * conn_values)         # [draw x group]
t_multi <- rep(NA_real_, nrow(theta_0))
