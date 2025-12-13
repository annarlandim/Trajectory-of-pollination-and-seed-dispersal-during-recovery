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
  
  

## important!! As long as I understand, this model assumes that theta can only have positive values because 
# the model is set for theta_inf and theta_old having positive values, which is the original scale of the Y_ old and Y_rec.
# in the case of composition, I will have negative values, so I need to think about that.
