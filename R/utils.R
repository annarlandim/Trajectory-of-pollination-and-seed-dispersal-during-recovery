hulls_hbt <- function(scores){
  
  hulls <- scores %>%
    group_by(Treatment3) %>%
    summarise(chull_indices = list(chull(RC1, RC2)))
  hull_points <- list()
  for (i in 1:nrow(hulls)) {
    indices <- c(hulls$chull_indices[[i]], hulls$chull_indices[[i]][1]) 
    points <- scores[scores$Treatment3 == hulls$Treatment3[i], ][indices, ]
    points$Treatment3 <- hulls$Treatment3[i] 
    hull_points[[i]] <- points
  }
  hull_points_df <- do.call(rbind, hull_points)
  return(hull_points_df)
}

calc_dist <- function(point1, point2) {
  sqrt(sum((point1 - point2) ^ 2))
}

originality <- function(unique_plot, scores){
  
  orig <- unique_plot %>%
    left_join(distinct(scores, !!sym(names(scores)[1]), !!sym(names(scores)[2]), .keep_all = TRUE), 
              by = c(names(scores)[1], names(scores)[2]), relationship = "many-to-many")
  
  centroids <- scores %>%
    group_by(Treatment3) %>%
    summarise(RC1 = mean(RC1), RC2 = mean(RC2))  
  
  orig$orig <- mapply(function(x, y, trat) {
    center <- centroids[centroids$Treatment3 == trat, c("RC1", "RC2")]
    calc_dist(c(x, y), center)
  }, orig$RC1, orig$RC2, orig$Treatment3)
  
  return(orig %>% mutate(count.t3 = NULL, count.plants.t3 = NULL))
}


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

recovery_tmulti_nimble <- function(samples, groups, conn, maxt = 1000) {
  # samples: mcmc.list (nimble output)
  # groups: integer indices dos grupos (ex: c(4,5,6))
  # conn: valor de conectividade (já padronizado), escalar
  # maxt: tempo máximo a procurar
  
  # junta chains em uma matriz: [draw x param]
  all_draws <- do.call(rbind, samples)
  cn <- colnames(all_draws)
  
  # nomes das colunas para os grupos selecionados
  theta0_cols   <- paste0("theta_0[",   groups, "]")
  thetainf_cols <- paste0("theta_inf[", groups, "]")
  alpha_cols    <- paste0("alpha_con[", groups, "]")
  beta_cols     <- paste0("beta_con[",  groups, "]")
  
  # extrai matrizes [draw x group]
  theta_0   <- as.matrix(all_draws[, theta0_cols,   drop = FALSE])
  theta_inf <- as.matrix(all_draws[, thetainf_cols, drop = FALSE])
  aMat      <- as.matrix(all_draws[, alpha_cols,    drop = FALSE])
  bMat      <- as.matrix(all_draws[, beta_cols,     drop = FALSE])
  
  # grade de tempo
  tgrid  <- seq(0, maxt, 0.1)
  
  # lambda por draw e grupo
  lambda <- exp(aMat + bMat * conn)   # [draw x group]
  
  # vetor de tempos de recuperação (um por draw)
  t_multi <- rep(NA_real_, nrow(theta_0))
  
  for (tt in tgrid) {
    # valor da métrica em cada tempo, draw e grupo
    At <- theta_0 + (theta_inf - theta_0) * (1 - exp(-lambda * tt))  # [draw x group]
    
    # desvio médio relativo em relação às assíntotas
    Mt <- rowMeans(abs(At - theta_inf) / theta_inf)
    
    hit <- is.na(t_multi) & (Mt <= 0.10)
    if (any(hit)) t_multi[hit] <- tt
    if (all(!is.na(t_multi))) break
  }
  
  t_multi
}

## weighted:

recovery_tmulti_nimble_w <- function(samples, groups, conn, weights, maxt = 1000) {
  all_draws <- do.call(rbind, samples)
  
  theta0_cols   <- paste0("theta_0[",   groups, "]")
  thetainf_cols <- paste0("theta_inf[", groups, "]")
  alpha_cols    <- paste0("alpha_con[", groups, "]")
  beta_cols     <- paste0("beta_con[",  groups, "]")
  
  theta_0   <- as.matrix(all_draws[, theta0_cols,   drop = FALSE])
  theta_inf <- as.matrix(all_draws[, thetainf_cols, drop = FALSE])
  aMat      <- as.matrix(all_draws[, alpha_cols,    drop = FALSE])
  bMat      <- as.matrix(all_draws[, beta_cols,     drop = FALSE])
  
  w <- as.numeric(weights)
  w <- w / sum(w)
  
  tgrid  <- seq(0, maxt, 0.1)
  lambda <- exp(aMat + bMat * conn)
  
  t_multi <- rep(NA_real_, nrow(theta_0))
  
  for (tt in tgrid) {
    At <- theta_0 + (theta_inf - theta_0) * (1 - exp(-lambda * tt))
    dev_rel <- abs(At - theta_inf) / theta_inf
    Mt <- as.numeric(dev_rel %*% w)
    
    hit <- is.na(t_multi) & (Mt <= 0.10)
    if (any(hit)) t_multi[hit] <- tt
    if (all(!is.na(t_multi))) break
  }
  
  t_multi
}

## without connectivity:

recovery_tmulti_nimble_w <- function(samples, groups, weights, maxt = 1000) {
  all_draws <- do.call(rbind, samples)
  
  theta0_cols   <- paste0("theta_0[",   groups, "]")
  thetainf_cols <- paste0("theta_inf[", groups, "]")
  lambda_cols   <- paste0("lambda[", groups, "]")
  
  theta_0   <- as.matrix(all_draws[, theta0_cols,   drop = FALSE])
  theta_inf <- as.matrix(all_draws[, thetainf_cols, drop = FALSE])
  lMat      <- as.matrix(all_draws[, lambda_cols,    drop = FALSE])
  
  w <- as.numeric(weights)
  w <- w / sum(w)
  
  tgrid  <- seq(0, maxt, 0.1)
  lambda <- lMat
  
  t_multi <- rep(NA_real_, nrow(theta_0))
  
  for (tt in tgrid) {
    At <- theta_0 + (theta_inf - theta_0) * (1 - exp(-lambda * tt))
    dev_rel <- abs(At - theta_inf) / theta_inf
    Mt <- as.numeric(dev_rel %*% w)
    
    hit <- is.na(t_multi) & (Mt <= 0.10)
    if (any(hit)) t_multi[hit] <- tt
    if (all(!is.na(t_multi))) break
  }
  
  t_multi
}

### per group:

recovery_tmulti_per_group_nimble <- function(samples, groups, conn,
                                             maxt = 1000) {
  
  # junta as chains
  all <- do.call(rbind, samples)
  
  # extrai draws [draw x group]
  theta_0   <- as.matrix(all[, paste0("theta_0[",   groups, "]"), drop = FALSE])
  theta_inf <- as.matrix(all[, paste0("theta_inf[", groups, "]"), drop = FALSE])
  aMat      <- as.matrix(all[, paste0("alpha_con[", groups, "]"), drop = FALSE])
  bMat      <- as.matrix(all[, paste0("beta_con[",  groups, "]"), drop = FALSE])
  
  n_draw  <- nrow(theta_0)
  n_group <- ncol(theta_0)
  
  lambda <- exp(aMat + bMat * conn)
  tgrid  <- seq(0, maxt, 0.1)
  
  # t_multi por draw por grupo
  tmat <- matrix(NA_real_, nrow = n_draw, ncol = n_group)
  colnames(tmat) <- paste0("g", groups)
  
  for (g in seq_len(n_group)) {
    T0   <- theta_0[, g]
    Tinf <- theta_inf[, g]
    lam  <- lambda[, g]
    
    done <- rep(FALSE, n_draw)
    
    for (tt in tgrid) {
      At  <- T0 + (Tinf - T0) * (1 - exp(-lam * tt))
      rel <- abs(At - Tinf) / Tinf
      
      hit <- (!done) & (rel <= 0.10)
      if (any(hit)) {
        tmat[hit, g] <- tt
        done[hit] <- TRUE
      }
      if (all(done)) break
    }
  }
  
  tmat
}

## without con

recovery_tmulti_per_group_nimble <- function(samples, groups, 
                                             maxt = 1000) {
  
  # junta as chains
  all <- do.call(rbind, samples)
  
  # extrai draws [draw x group]
  theta_0   <- as.matrix(all[, paste0("theta_0[",   groups, "]"), drop = FALSE])
  theta_inf <- as.matrix(all[, paste0("theta_inf[", groups, "]"), drop = FALSE])
  lMat      <- as.matrix(all[, paste0("lambda[", groups, "]"), drop = FALSE])
  
  n_draw  <- nrow(theta_0)
  n_group <- ncol(theta_0)
  
  lambda <- lMat
  tgrid  <- seq(0, maxt, 0.1)
  
  # t_multi por draw por grupo
  tmat <- matrix(NA_real_, nrow = n_draw, ncol = n_group)
  colnames(tmat) <- paste0("g", groups)
  
  for (g in seq_len(n_group)) {
    T0   <- theta_0[, g]
    Tinf <- theta_inf[, g]
    lam  <- lambda[, g]
    
    done <- rep(FALSE, n_draw)
    
    for (tt in tgrid) {
      At  <- T0 + (Tinf - T0) * (1 - exp(-lam * tt))
      rel <- abs(At - Tinf) / Tinf
      
      hit <- (!done) & (rel <= 0.10)
      if (any(hit)) {
        tmat[hit, g] <- tt
        done[hit] <- TRUE
      }
      if (all(done)) break
    }
  }
  
  tmat
}

