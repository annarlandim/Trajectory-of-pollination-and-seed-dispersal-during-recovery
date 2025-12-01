rec_conn <- dataSub$ConIndex[dataSub$type == "rec"]
conn_values <- c(
  Low    = quantile(rec_conn, 0.25, na.rm = TRUE),
  Medium = quantile(rec_conn, 0.50, na.rm = TRUE),
  High   = quantile(rec_conn, 0.75, na.rm = TRUE)
)


tmultiList_conn <- list(
  SD = lapply(conn_values, function(cn) recovery_tmulti(samp1, groups = 4:6, conn = cn)),
  FS = lapply(conn_values, function(cn) recovery_tmulti(samp1, groups = 1:3, conn = cn))
)

theta_inf <- do.call(rbind, as.mcmc.list(samp1$theta_inf))
theta_0 <- do.call(rbind, as.mcmc.list(samp1$theta_0))
alpha_con <- do.call(rbind, as.mcmc.list(samp1$alpha_con))
beta_con <- do.call(rbind, as.mcmc.list(samp1$beta_con))

# sigma_old <- do.call(rbind, as.mcmc.list(samp1$sigma_old))
# sigma_rec <- do.call(rbind, as.mcmc.list(samp1$sigma_rec))
# sigmaSq_old <- sigma_old^2
# sigmaSq_rec <- sigma_rec^2

sigmaSq_old <- do.call(rbind, as.mcmc.list(samp1$sigmaSq_old))
sigmaSq_rec <- do.call(rbind, as.mcmc.list(samp1$sigmaSq_rec))

summ <- function(x) quantile(as.numeric(x), probs = c(.05,.25,.5,.75,.95), na.rm = TRUE)

qtMulti_conn <- lapply(tmultiList_conn, function(proc_list) {
  lapply(proc_list, summ)
})

t_multi_mat <- cbind(
  SD_Low    = as.numeric(tmultiList_conn$SD$Low),
  SD_Medium = as.numeric(tmultiList_conn$SD$Medium),
  SD_High   = as.numeric(tmultiList_conn$SD$High),
  FS_Low    = as.numeric(tmultiList_conn$FS$Low),
  FS_Medium = as.numeric(tmultiList_conn$FS$Medium),
  FS_High   = as.numeric(tmultiList_conn$FS$High)
)

qtheta_0   <- apply(theta_0,   2, quantile, prob = c(.05,.25,.5,.75,.95))
qtheta_inf <- apply(theta_inf, 2, quantile, prob = c(.05,.25,.5,.75,.95))
qalpha     <- apply(alpha_con,     2, quantile, prob = c(.05,.25,.5,.75,.95))
qbeta      <- apply(beta_con,      2, quantile, prob = c(.05,.25,.5,.75,.95))
qsigmarec  <- apply(sigmaSq_rec, 2, quantile, prob = c(.05,.25,.5,.75,.95))
qsigmaold  <- apply(sigmaSq_old, 2, quantile, prob = c(.05,.25,.5,.75,.95))

#### per group ####

# Per-group t90 for selected groups at a given connectivity
# groups: integer indices of columns in Y (e.g., 4:6 for bats/birds/nf mammals)
# conn:   scalar connectivity
# maxt/by: time grid

recovery_t90_per_group <- function(samples, groups, conn, maxt = 1000) {
  theta_0   <- do.call(rbind, as.mcmc.list(samples$theta_0))[ , paste0("theta_0[",   groups, "]"), drop = FALSE]
  theta_inf   <- do.call(rbind, as.mcmc.list(samples$theta_inf))[ , paste0("theta_inf[",   groups, "]"), drop = FALSE]
  aMat   <- do.call(rbind, as.mcmc.list(samples$alpha_con))[ , paste0("alpha_con[",   groups, "]"), drop = FALSE]
  bMat   <- do.call(rbind, as.mcmc.list(samples$beta_con))[ , paste0("beta_con[",   groups, "]"), drop = FALSE]
  
  n_draw  <- nrow(theta_0)
  n_group <- ncol(theta_0)
  tgrid   <- seq(0, maxt, by = 0.1)
  
  # λ per draw & group
  lambda <- exp(aMat + bMat * conn)   # [draw x group]
  
  # result matrix: one t90 per draw per group
  t90 <- matrix(NA_real_, nrow = n_draw, ncol = n_group,
                dimnames = list(NULL, paste0("g", groups)))
  
  # loop groups (fast, keeps memory modest)
  for (g in seq_len(n_group)) {
    T0   <- theta_0[ , g]
    Tinf <- theta_inf[, g]
    lam  <- lambda[, g]
    
    # Track which draws already met criterion
    filled <- rep(FALSE, n_draw)
    
    for (tt in tgrid) {
      At <- T0 + (Tinf - T0) * (1 - exp(-lam * tt))
      # increasing vs decreasing criterion
      inc <- (T0 <= Tinf)
      hit_inc <-  inc & !filled & (At >= 0.9 * Tinf)
      hit_dec <- !inc & !filled & (At <= 1.1 * Tinf)
      hit <- hit_inc | hit_dec
      if (any(hit)) {
        t90[hit, g] <- tt
        filled[hit] <- TRUE
      }
      if (all(filled)) break
    }
  }
  t90
}

rec_conn <- dataSub$ConIndex[dataSub$type == "rec"]

c_low  <- quantile(rec_conn, 0.25, na.rm = TRUE)
c_mid  <- quantile(rec_conn, 0.50, na.rm = TRUE)
c_high <- quantile(rec_conn, 0.75, na.rm = TRUE)

t90_SD_low  <- recovery_t90_per_group(samp1, groups = 4:6, conn = c_low)
t90_SD_mid  <- recovery_t90_per_group(samp1, groups = 4:6, conn = c_mid)
t90_SD_high <- recovery_t90_per_group(samp1, groups = 4:6, conn = c_high)

t90_FS_low  <- recovery_t90_per_group(samp1, groups = 1:3, conn = c_low)
t90_FS_mid  <- recovery_t90_per_group(samp1, groups = 1:3, conn = c_mid)
t90_FS_high <- recovery_t90_per_group(samp1, groups = 1:3, conn = c_high)

summ <- function(x) quantile(as.numeric(x), probs = c(.05,.25,.5,.75,.95), na.rm = TRUE)

apply(t90_SD_low,  2, summ)   # Low connectivity
apply(t90_SD_mid,  2, summ)   # Median
apply(t90_SD_high, 2, summ)   # High

apply(t90_FS_low,  2, summ)   # Low connectivity
apply(t90_FS_mid,  2, summ)   # Median
apply(t90_FS_high, 2, summ)   # High
