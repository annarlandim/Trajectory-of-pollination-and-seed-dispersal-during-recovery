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


hulls_groups <- function(scores){
  
  hulls <- scores %>%
    group_by(group) %>%
    summarise(chull_indices = list(chull(RC1, RC2)))
  hull_points <- list()
  for (i in 1:nrow(hulls)) {
    indices <- c(hulls$chull_indices[[i]], hulls$chull_indices[[i]][1]) 
    points <- scores[scores$group == hulls$group[i], ][indices, ]
    points$group <- hulls$group[i] 
    hull_points[[i]] <- points
  }
  hull_points_df <- do.call(rbind, hull_points)
  return(hull_points_df)
}

hulls_groups_stage <- function(scores){
  hulls <- scores %>%
    group_by(group, Treatment3) %>%
    summarise(chull_indices = list(chull(RC1, RC2)), .groups = "drop")
  hull_points <- list()
  for (i in 1:nrow(hulls)) {
    indices <- c(hulls$chull_indices[[i]], hulls$chull_indices[[i]][1])
    subset_pts <- scores[scores$group == hulls$group[i] & scores$Treatment3 == hulls$Treatment3[i], ]
    points <- subset_pts[indices, ]
    points$group <- hulls$group[i]
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

gelman_log <- function(samples_list, pattern) {
  
  names_vars <- varnames(samples_list)
  idx <- grep(paste0(pattern, "\\["), names_vars) 
  
  samples_log <- as.mcmc.list(lapply(samples_list, function(chain) {
    
    subset_chain <- chain[, idx, drop = FALSE]
    
    return(as.mcmc(log(subset_chain))) 
  }))
  return(gelman.diag(samples_log))
}

compute_multi_trajectory <- function(draw, idxs, weights, conn, t) {
  
  At <- numeric(length(idxs))
  
  for(i in seq_along(idxs)) {
    j <- idxs[i]
    
    t0   <- draw[paste0("theta_0[", j, "]")]
    tinf <- draw[paste0("theta_inf[", j, "]")]
    a    <- draw[paste0("alpha_con[", j, "]")]
    b    <- draw[paste0("beta_con[", j, "]")]
    
    lambda <- exp(a + b * conn)
    
    At[i] <- t0 + (tinf - t0) * (1 - exp(-lambda * t))
  }
  
  sum(At * weights)
}

plot_trajectory_multi <- function(var_indices, weights,
                                  conn, title, colors,
                                  data_subset, leg.pos,
                                  names_groups) {
  
  plot(NULL, xlim = c(0, 40), ylim = c(0.2, 0.8), 
       xlab = "Time (years)", ylab = "Functional Diversity",
       bty = "l", cex.lab = 1.2)
  mtext(title, side = 3, line = 1.5, adj = 0, font = 2, cex = 1)
  
  # ----- 1) Individual groups  -----
  
  for(i in seq_along(var_indices)) {
    idx <- var_indices[i]
    
    sub_lines <- sample(1:nrow(mcmc_mat), 50)
    
    for(s in sub_lines) {
      t0 <- mcmc_mat[s, paste0("theta_0[", idx, "]")]
      tinf <- mcmc_mat[s, paste0("theta_inf[", idx, "]")]
      a <- mcmc_mat[s, paste0("alpha_con[", idx, "]")]
      b <- mcmc_mat[s, paste0("beta_con[", idx, "]")]
      
      lambda <- exp(a + b * conn)
      
      curve(t0 + (tinf - t0) * (1 - exp(-lambda * x)), 
            add = TRUE,
            col = adjustcolor(colors[i], alpha.f = 0.3),
            lwd = 0.5)
    }
    
    # Median line
    t0_m   <- median(mcmc_mat[, paste0("theta_0[", idx, "]")])
    tinf_m <- median(mcmc_mat[, paste0("theta_inf[", idx, "]")])
    a_m    <- median(mcmc_mat[, paste0("alpha_con[", idx, "]")])
    b_m    <- median(mcmc_mat[, paste0("beta_con[", idx, "]")])
    
    lamb_m <- exp(a_m + b_m * conn)
    
    curve(t0_m + (tinf_m - t0_m) * (1 - exp(-lamb_m * x)), 
          add = TRUE, col = colors[i], lwd = 3)
    
  }
  
  # ----- 2) Across groups -----
  
  t_seq <- seq(0, 40, length = 200)
  
  # For each time and draw
  multi_mat <- sapply(t_seq, function(tt) {
    apply(mcmc_mat, 1, function(draw)
      compute_multi_trajectory(draw, var_indices, weights, conn, tt)
    )
  })
  
  # Summaries
  med  <- apply(multi_mat, 2, median)
  loqt <- apply(multi_mat, 2, quantile, 0.025)
  hiqt <- apply(multi_mat, 2, quantile, 0.975)
  
  # Ribbon
  polygon(c(t_seq, rev(t_seq)),
          c(loqt, rev(hiqt)),
          col = adjustcolor(colors[4], 0.15),
          border = NA)
  
  # Median line 
  lines(t_seq, med, col = colors[4], lwd = 4)
  
  legend(leg.pos, title = "Interactions driven by:",
         legend = c(names_groups, "Across groups"),
         col = c(colors, "black"),
         lwd = c(rep(3, length(colors)), 4),
         bty = "n", y.intersp = 0.7)
}

plot_old_growth_multi <- function(var_indices, weights, colors) {
  
  n <- length(var_indices)
  
  plot(NULL,
       xlim = c(0.5, n + 1.5),
       ylim = c(0.2, 0.8), 
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = "",
       bty = "n")
  
  axis(1,
       at = c(0.5, n + 1.5),
       labels = FALSE,
       lwd.ticks = 0, lwd = 1)
  
  # ----- 1) Individual groups -----
  
  for(i in seq_along(var_indices)) {
    
    idx <- var_indices[i]
    
    vals <- mcmc_mat[, paste0("theta_inf[", idx, "]")]
    
    stats <- quantile(vals,
                      probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    
    segments(i, stats[1], i, stats[5],
             col = colors[i], lwd = 1)
    
    segments(i, stats[2], i, stats[4],
             col = colors[i], lwd = 4)
    
    points(i, stats[3],
           pch = 21, bg = colors[i],
           col = colors[i],
           cex = 1.5, lwd = 2)
  }
  
  # ----- 2) MULTIFUNCTIONALITY -----
  
  # Compute posterior draws of multi asymptote
  multi_vals <- apply(mcmc_mat, 1, function(draw) {
    
    thetas <- sapply(var_indices, function(j)
      draw[paste0("theta_inf[", j, "]")]
    )
    
    sum(thetas * weights)
  })
  
  stats_m <- quantile(multi_vals,
                      probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  
  i_multi <- n + 1
  
  segments(i_multi, stats_m[1], i_multi, stats_m[5],
           col = colors[4], lwd = 1)
  
  segments(i_multi, stats_m[2], i_multi, stats_m[4],
           col = colors[4], lwd = 4)
  
  points(i_multi, stats_m[3],
         pch = 21, bg = colors[4],
         col = colors[4],
         cex = 1.7, lwd = 2)
  
  
  axis(1,
       at = 1:(n+1),
       labels = FALSE,
       lwd.ticks = 0,
       las = 2, cex.axis = 0.8)
  
  mtext("Old-growth", side = 1, line = 2.1, cex = 0.7)
}

get_net_change <- function(j) {
  theta_0   <- mcmc_mat[, paste0("theta_0[",   j, "]")]
  theta_inf <- mcmc_mat[, paste0("theta_inf[", j, "]")]
  theta_inf - theta_0
}

get_multi_net_change <- function(idxs, weights) {
  out <- numeric(nrow(mcmc_mat))
  for (i in seq_along(idxs)) {
    out <- out + weights[i] * get_net_change(idxs[i])
  }
  out
}

ci_summary <- function(x) quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)

draw_netchange_plot <- function(current_list, label, x_lim) {
  
  total_height <- sum(sapply(current_list, function(x) if (!is.null(x$height)) x$height else 1))
  
  plot(NULL, xlim = x_lim, ylim = c(0.5, total_height + 0.5),
       xlab = "Net change", ylab = "",
       yaxt = "n", bty = "l")
  
  mtext(label, side = 3, line = 1.5, adj = 0, font = 2, cex = 1.2)
  
  abline(v = 0, lty = 2, col = "grey50")   # no net change reference
  
  y_tracker <- total_height + 0.5
  
  for (i in seq_along(current_list)) {
    item  <- current_list[[i]]
    h     <- if (!is.null(item$height)) item$height else 1
    y_pos <- y_tracker - (h / 2)
    y_tracker <- y_tracker - h
    
    if (is.null(item$type) || item$type == "spacer") next
    
    stats <- item$val
    lwd_thin  <- ifelse(item$type == "multi", 2, 2)
    lwd_thick <- ifelse(item$type == "multi", 6, 6)
    
    segments(stats[1], y_pos, stats[5], y_pos, col = item$col, lwd = lwd_thin)   # 90% CI
    segments(stats[2], y_pos, stats[4], y_pos, col = item$col, lwd = lwd_thick)  # IQR
    points(stats[3], y_pos, pch = 21, bg = item$col, col = "white",
           cex = ifelse(item$type == "multi", 2.8, 2.5))
    
    axis(2, at = y_pos, labels = item$lab, las = 1, cex.axis = 0.9,
         font = ifelse(item$type == "multi", 2, 1), tick = TRUE, tck = -0.02)
  }
}

recovery_tmulti_nimble_w_con <- function(samples, groups, conn, weights, maxt = 1000) {
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

### per group:

recovery_tmulti_per_group_nimble_con <- function(samples, groups, conn,
                                             maxt = 1000) {
  
  all <- do.call(rbind, samples)
  
  theta_0   <- as.matrix(all[, paste0("theta_0[",   groups, "]"), drop = FALSE])
  theta_inf <- as.matrix(all[, paste0("theta_inf[", groups, "]"), drop = FALSE])
  aMat      <- as.matrix(all[, paste0("alpha_con[", groups, "]"), drop = FALSE])
  bMat      <- as.matrix(all[, paste0("beta_con[",  groups, "]"), drop = FALSE])
  
  n_draw  <- nrow(theta_0)
  n_group <- ncol(theta_0)
  
  lambda <- exp(aMat + bMat * conn)
  tgrid  <- seq(0, maxt, 0.1)
  
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

diagnose_block <- function(block, samples, samples_all) {
  
  names_vars <- varnames(samples)
  idx <- grep(block$pattern, names_vars)
  
  if (length(idx) == 0) {
    warning(sprintf("No parameters matched for block '%s'", block$name))
    return(NULL)
  }
  
  if (block$log_scale) {
    # Gelman on log scale 
    samples_list <- as.mcmc.list(lapply(samples, function(chain) {
      as.mcmc(log(chain[, idx, drop = FALSE]))
    }))
    ess_input <- log(samples_all[, idx, drop = FALSE])
  } else {
    samples_list <- samples[, idx, drop = FALSE]
    ess_input <- samples_all[, idx, drop = FALSE]
  }
  
  gd <- gelman.diag(samples_list, autoburnin = FALSE, multivariate = FALSE)
  rhat_vals <- gd$psrf[, "Point est."]
  
  ess_vals <- effectiveSize(ess_input)
  
  tibble(
    block = block$name,
    parameter = names(rhat_vals),
    rhat = as.numeric(rhat_vals),
    ess = as.numeric(ess_vals[names(rhat_vals)])
  )
}

summarise_param <- function(base_name, pretty_name, samples, samples_all,
                            indices = 1:6, transform_fn = identity) {
  
  pattern <- paste0("^", base_name, "\\[")
  names_vars <- varnames(samples)
  idx_cols <- grep(pattern, names_vars)
  
  col_idx_num <- as.integer(str_extract(names_vars[idx_cols], "(?<=\\[)\\d+(?=\\])"))
  keep <- idx_cols[col_idx_num %in% indices]
  keep <- keep[order(col_idx_num[col_idx_num %in% indices])]
  
  draws <- samples_all[, keep, drop = FALSE]
  draws_t <- apply(draws, 2, transform_fn)
  
  q <- t(apply(draws_t, 2, quantile, probs = c(0.5, 0.025, 0.25, 0.75, 0.975)))
  colnames(q) <- c("Median", "q2.5", "q25", "q75", "q97.5")
  
  neff <- effectiveSize(as.mcmc(draws_t))
  
  samples_t_list <- as.mcmc.list(lapply(samples, function(chain) {
    sub <- chain[, keep, drop = FALSE]
    as.mcmc(apply(sub, 2, transform_fn))
  }))
  psrf <- gelman.diag(samples_t_list, autoburnin = FALSE, multivariate = FALSE)$psrf[, "Point est."]
  
  tibble(
    Parameter = pretty_name,
    Group = group_labels[as.character(sort(indices))],
    Median = q[, "Median"], `2.5%` = q[, "q2.5"], `25%` = q[, "q25"],
    `75%` = q[, "q75"], `97.5%` = q[, "q97.5"],
    Neff = as.numeric(neff), PSRF = as.numeric(psrf)
  )
}

tau <- function(sigma) 1 / sigma^2