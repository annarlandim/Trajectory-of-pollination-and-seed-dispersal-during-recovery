#### Data
## With connectivity
qt90 <- readRDS("output/t90_con.rds")
samples <- readRDS("output/model_posteriors_con.rds")
mcmc_mat <- as.matrix(samples)

weights_df <- read.csv("data/processed/weights_df.csv")

data <- read.csv("data/processed/model_df.csv")
data$type <- factor(ifelse(1:nrow(data) %in% grep("OG", data$Plot_ID), "old", "rec"),
                    levels = c("old", "rec"))


data <- data %>% mutate(ConIndex = as.numeric(scale(ConIndex)))

rec_conn <- data$ConIndex[data$type == "rec"]
conn_values <- c(
  Low    = quantile(rec_conn, 0.25, na.rm = TRUE),
  Medium = quantile(rec_conn, 0.50, na.rm = TRUE),
  High   = quantile(rec_conn, 0.75, na.rm = TRUE)
)

dataSub <- subset(data, select = c(type, RegTime,  
                                   FDBees, FDMoths,  FDBat_pol, 
                                   FDBats, FDBirds, FDNf#, 
                                   # FDSeeds, AbSeeds, RichSeeds, ShnSeeds,
                                   # FDSdlng, AbSdlng, RichSdlng, ShnSdlng
                                   )) 
dataSub <- dataSub %>%
  mutate(across(
    .cols = -c(1:2),     
    # scaling by dividing the max per column.
    # this is the best for logscale, where negatives are not allowed
    .fns = ~ .x / max(.x, na.rm = TRUE) 
  ))

str(dataSub)

zeros <- dataSub %>% 
  map_lgl(~ any(. == 0, na.rm = TRUE))
print(zeros)

# to avoid problems with log(0)
# eps <- 1e-6
# dataSub$AbSdlng <- pmax(dataSub$AbSdlng, eps)
# dataSub$RichSdlng <- pmax(dataSub$RichSdlng, eps)
# dataSub$ShnSdlng <- pmax(dataSub$ShnSdlng, eps)

long <- dataSub %>%
  pivot_longer(cols = -c(1:2), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = colnames(dataSub)[-c(1:2)]), variable = as.integer(variable)) %>%
  filter(!is.na(value))

### Recovery trajectories

var_names <- colnames(dataSub)[-c(1,2)]

#### Multifunctionality

cols_pol <- c("#c994c7", "#9e9ac8", "#cbc9e2", "#3f007d")
cols_sd <- c("#3182bd", "#7bccc4", "#bdd7e7", "#08306b")


# Calculate Weighted Data Points:

w_pol <- c(FDBees=weights_df[1,"prctg"]/sum(weights_df[1:3,"prctg"]), FDMoths=weights_df[2,"prctg"]/sum(weights_df[1:3,"prctg"]), FDBat_pol=weights_df[3,"prctg"]/sum(weights_df[1:3,"prctg"])) 
w_sd <- c(FDBats=weights_df[4,"prctg"]/sum(weights_df[4:6,"prctg"]), FDBirds=weights_df[5,"prctg"]/sum(weights_df[4:6,"prctg"]), FDNf=weights_df[6,"prctg"]/sum(weights_df[4:6,"prctg"]))

dataSub$multi_pol <- apply(dataSub[, names(w_pol)], 1, function(x) {
  sum(x * w_pol, na.rm = TRUE)
})
dataSub$multi_sd <-apply(dataSub[, names(w_sd)], 1, function(x) {
  sum(x * w_sd, na.rm = TRUE)
})

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
  
  # ----- 1) Individual groups (as you already had) -----
  
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
    
    # points(data_subset$RegTime[data_subset$variable == idx], 
    #        data_subset$value[data_subset$variable == idx], 
    #        pch = 21,
    #        bg = adjustcolor(colors[i], alpha.f = 1),
    #        col = "white", cex = 1.2)
  }
  
  # ----- 2) MULTIFUNCTIONALITY -----
  
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
  
  # Median line multifunctionality
  lines(t_seq, med, col = colors[4], lwd = 4)
  
  legend(leg.pos, title = "Interactions driven by:",
         legend = c(names_groups, "Multifunctionality"),
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

# reordering for plotting

w_pol <- w_pol[c("FDBat_pol", "FDBees", "FDMoths")]

# Attention: Pollinators are ordered as 3, 1, 2
names_pol <- c("Bats", "Bees", "Moths")
names_sd <- c("Bats", "Birds", "Non-flying mammals")

# svg("output/Figures/Figure_trends_ts_unique.7.svg", width = 8, height = 5)

layout_mat <- matrix(c(1, 2, 3, 4), nrow = 1)
layout(layout_mat, widths = c(3, 0.8, 3, 0.8))

par(mar = c(4, 4, 3, 0.8), oma = c(1, 1, 1, 1), mgp = c(2, 0.7, 0), tcl = -0.3)

# --- PAINEL A: POLLINATION ---
par(mar = c(4, 4, 3, 0.8))
plot_trajectory_multi(c(3,1,2), w_pol, conn_values["Medium.50%"], "A) Pollination", cols_pol, long[long$type == "rec", ], "bottomleft", names_pol)
par(mar = c(4, 0.8, 3, 1))
plot_old_growth_multi(c(3,1,2), w_pol, cols_pol)

# --- PAINEL B: SEED DISPERSAL ---
par(mar = c(4, 4, 3, 0.8))
plot_trajectory_multi(4:6, w_sd, conn_values["Medium.50%"], "B) Seed dispersal", cols_sd, long[long$type == "rec", ], "bottomright", names_sd)
par(mar = c(4, 0.8, 3, 1))
plot_old_growth_multi(4:6, w_sd, cols_sd)

# dev.off()

# Figure 3: Recovery time 

# --- Pollination List (Panel A) ---
pollination_list <- list(
  list(val = qt90$pollination$per_group$`Medium.50%`[,3], lab = "Bats", col = cols_pol[1], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = qt90$pollination$per_group$`Medium.50%`[,1], lab = "Bees", col = cols_pol[2], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = qt90$pollination$per_group$`Medium.50%`[,2], lab = "Moths", col = cols_pol[3], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = qt90$pollination$multi$`Medium.50%`, lab = "Multi\nfunctionality", col = cols_pol[4], type ="multi")
)

# --- Seed Dispersal List (Panel B) ---
dispersal_list <- list(
  list(val = qt90$dispersal$per_group$`Medium.50%`[,1], lab = "Bats", col = cols_sd[1], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = qt90$dispersal$per_group$`Medium.50%`[,2], lab = "Birds", col = cols_sd[2], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = qt90$dispersal$per_group$`Medium.50%`[,3], lab = "Non-flying\nmammals", col = cols_sd[3], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = qt90$dispersal$multi$`Medium.50%`, lab = "Multi\nfunctionality", col = cols_sd[4], type ="multi")
)

draw_rectime_plot <- function(current_list, label, x_max_val) {
  total_height <- sum(sapply(current_list, function(x) if(!is.null(x$height)) x$height else 1))
  
  plot(NULL, xlim = c(0, x_max_val), ylim = c(0.5, total_height + 0.5), 
       xlab = "Time to convergence (years)", ylab = "", 
       yaxt = "n", bty = "l")
  
  mtext(label, side = 3, line = 1.5, adj = 0, font = 2, cex = 1.2)
  
  y_tracker <- total_height + 0.5
  
  for(i in seq_along(current_list)) {
    item <- current_list[[i]]
    h <- if(!is.null(item$height)) item$height else 1
    y_pos <- y_tracker - (h / 2)
    y_tracker <- y_tracker - h
    
    if(is.null(item$type) || item$type == "spacer") next 
    
    stats <- item$val
    # Error Bars
    segments(stats[1], y_pos, stats[5], y_pos, col = item$col, lwd = 2) # 90% CI
    segments(stats[2], y_pos, stats[4], y_pos, col = item$col, lwd = 6) # IQR
    points(stats[3], y_pos, pch = 21, bg = item$col, col = "white", cex = 2.5)
    
    # Y Axis Labels
    axis(2, at = y_pos, labels = item$lab, las = 1, cex.axis = 0.9, 
         font = 1, tick = TRUE, tck = -0.02)
  }
}

# Use a consistent x-axis across both for easy comparison
pol_x_max <- max(as.numeric(qt90$pol$per_group$`Medium.50%`)) * 1.1

sd_x_max <- max(as.numeric(qt90$dispersal$per_group$`Medium.50%`)) * 1.1

# svg("output/Figures/Figure_rectime_ts_unique.svg", width = 10, height = 5)
# Set up layout for two panels (one above the other or side-by-side)
# par(mfrow = c(1, 2), mar = c(5, 8, 3, 2), oma = c(0, 0, 0, 0))

layout(matrix(c(1, 2), nrow = 1), widths = c(1+2, 4+2))

# Adjust margins: You may want to decrease the left margin for Panel B 
# since it no longer needs as much space for labels if they are short.
par(mar = c(5, 6, 3, 1))
# Draw Panel A
draw_rectime_plot(pollination_list, "A) Pollination", pol_x_max)

# Draw Panel B
par(mar = c(5, 4, 3, 2))
draw_rectime_plot(dispersal_list, "B) Seed Dispersal", sd_x_max)

# dev.off()

################# old stuff ############################################
## New figure 3:

get_group_recovery_traj <- function(samples, groups, conn, x_lim) {
  
  all_draws <- do.call(rbind, samples)
  
  res <- list()
  
  for(g in groups) {
    
    theta_0   <- all_draws[, paste0("theta_0[", g, "]")]
    theta_inf <- all_draws[, paste0("theta_inf[", g, "]")]
    a         <- all_draws[, paste0("alpha_con[", g, "]")]
    b         <- all_draws[, paste0("beta_con[", g, "]")]
    
    lambda <- exp(a + b * conn)
    
    traj <- matrix(NA, length(theta_0), length(x_lim))
    
    for(j in seq_along(x_lim)) {
      
      tt <- x_lim[j]
      
      At <- theta_0 + (theta_inf - theta_0) * (1 - exp(-lambda * tt))
      dev_rel <- abs(At - theta_inf) / theta_inf
      
      traj[, j] <- 1 - dev_rel
      
    }
    
    res[[paste0("g", g)]] <- traj
  }
  
  res
}

draw_recovery <- function(traj, col, x_lim) {
  
  # This doesnt make sense because is confidence, not credible interval
  # polygon(
  #   c(x_lim, rev(x_lim)),
  #   c(apply(traj, 2, quantile, 0.025),
  #     rev(apply(traj, 2, quantile, 0.975))),
  #   col = adjustcolor(col, 0.15),
  #   border = NA
  # )
  
  lines(x_lim, apply(traj, 2, median), col = col, lwd = 3)
}

plot_recovery_panel <- function(groups, x_lim, weights, conn, colors, title) {
  
  plot(NULL,
       xlim = range(x_lim),
       ylim = c(0,1),
       xlab = "Time (years)",
       ylab = "Relative convergence",
       bty = "l")
  
  mtext(title, side = 3, adj = 0, font = 2, cex = 1.2)
  
  # individual groups
  group_traj <- get_group_recovery_traj(samples, groups, conn, x_lim)
  
  for(i in seq_along(groups)){
    draw_recovery(group_traj[[i]], colors[i], x_lim)
  }
  
  # multifunctionality
  multi_traj <- get_multi_recovery_traj(
    samples,
    groups,
    conn,
    weights,
    x_lim 
  )
  
  draw_recovery(multi_traj, colors[4], x_lim)
  
  abline(h=0.9, lty=2, col="grey40")
}

svg("output/Figures/Figure_convergence.svg", width = 8, height = 5)

layout(matrix(c(1,2), nrow = 1))

par(mar = c(4,4,3,1),
    oma = c(1,1,1,1),
    mgp = c(2,0.7,0),
    tcl = -0.3)

# PANEL A
plot_recovery_panel(
  groups = c(3,1,2),
  x_lim = seq(0, 25, length = 200),
  weights = w_pol,
  conn = conn_values["Medium.50%"],
  colors = cols_pol,
  title = "A) Pollination"
)

# PANEL B
plot_recovery_panel(
  groups = 4:6,
  x_lim = seq(0, 110, length = 200),
  weights = w_sd,
  conn = conn_values["Medium.50%"],
  colors = cols_sd,
  title = "B) Seed dispersal"
)

dev.off()

## Both together:

svg("output/Figures/Figure3.svg", width = 10, height = 8)

layout(matrix(c(1,2,
                3,4), nrow = 2, byrow = TRUE))

par(oma = c(1,1,1,1),
    mgp = c(2,0.7,0),
    tcl = -0.3)

par(mar = c(4,4,3,1))

plot_recovery_panel(
  groups = c(3,1,2),
  x_lim = seq(0, 25, length = 200),
  weights = w_pol,
  conn = conn_values["Medium.50%"],
  colors = cols_pol,
  title = "A)"
)

par(mar = c(4,4,3,1))

plot_recovery_panel(
  groups = 4:6,
  x_lim = seq(0, 110, length = 200),
  weights = w_sd,
  conn = conn_values["Medium.50%"],
  colors = cols_sd,
  title = "B)"
)

par(mar = c(5,6,2,1))

draw_rectime_plot(
  pollination_list,
  "C)",
  pol_x_max
)

par(mar = c(5,4,2,2))

draw_rectime_plot(
  dispersal_list,
  "D)",
  sd_x_max
)

dev.off()

###########################################################################################

par(mar = c(4, 10, 3, 1))

# summary_list <- list(
#   
#   # Pollination
#   
#   #High
#   list(val = qt90$pollination$multi$`High.75%`, lab = "", col = cols_pol[4], type ="multi"),
#   list(val = qt90$pollination$per_group$`High.75%`[,1], lab = "", col = cols_pol[1], type = "group"),
#   list(val = qt90$pollination$per_group$`High.75%`[,2], lab = "", col = cols_pol[2], type = "group"),
#   list(val = qt90$pollination$per_group$`High.75%`[,3], lab = "", col = cols_pol[3], type = "group"),
#   list(type = "spacer", height = 0.5),
#   #Medium
#   list(val = qt90$pollination$multi$`Medium.50%`, lab = "", col = cols_pol[4], type ="multi"),
#   list(val = qt90$pollination$per_group$`Medium.50%`[,1], lab = "", col = cols_pol[1], type = "group"),
#   list(val = qt90$pollination$per_group$`Medium.50%`[,2], lab = "", col = cols_pol[2], type = "group"),
#   list(val = qt90$pollination$per_group$`Medium.50%`[,3], lab = "", col = cols_pol[3], type = "group"),
#   list(type = "spacer", height = 0.5),
#   #Low
#   list(val = qt90$pollination$multi$`Low.25%`, lab = "", col = cols_pol[4], type ="multi"),
#   list(val = qt90$pollination$per_group$`Low.25%`[,1], lab = "", col = cols_pol[1], type = "group"),
#   list(val = qt90$pollination$per_group$`Low.25%`[,2], lab = "", col = cols_pol[2], type = "group"),
#   list(val = qt90$pollination$per_group$`Low.25%`[,3], lab = "", col = cols_pol[3], type = "group"),
#   
#   # Visual space
#   list(type = "spacer", height = 2),
#   
#   # Seed dispersal
#   
#   #High
#   list(val = qt90$dispersal$multi$`High.75%`, lab = "", col = cols_sd[4], type ="multi"),
#   list(val = qt90$dispersal$per_group$`High.75%`[,1], lab = "", col = cols_sd[1], type = "group"),
#   list(val = qt90$dispersal$per_group$`High.75%`[,2], lab = "", col = cols_sd[2], type = "group"),
#   list(val = qt90$dispersal$per_group$`High.75%`[,3], lab = "", col = cols_sd[3], type = "group"),
#   list(type = "spacer", height = 0.5),
#   #Medium
#   list(val = qt90$dispersal$multi$`Medium.50%`, lab = "", col = cols_sd[4], type ="multi"),
#   list(val = qt90$dispersal$per_group$`Medium.50%`[,1], lab = "", col = cols_sd[1], type = "group"),
#   list(val = qt90$dispersal$per_group$`Medium.50%`[,2], lab = "", col = cols_sd[2], type = "group"),
#   list(val = qt90$dispersal$per_group$`Medium.50%`[,3], lab = "", col = cols_sd[3], type = "group"),
#   list(type = "spacer", height = 0.5),
#   #Low
#   list(val = qt90$dispersal$multi$`Low.25%`, lab = "", col = cols_sd[4], type ="group"),
#   list(val = qt90$dispersal$per_group$`Low.25%`[,1], lab = "", col = cols_sd[1], type = "group"),
#   list(val = qt90$dispersal$per_group$`Low.25%`[,2], lab = "", col = cols_sd[2], type = "group"),
#   list(val = qt90$dispersal$per_group$`Low.25%`[,3], lab = "", col = cols_sd[3], type = "group")
# )

summary_list <- list(
  
  # Pollination
  
  #Bees
  # list(val = qt90$pollination$per_group$`High.75%`[,1], lab = "H", col = cols_pol[1], type = "group"),
  list(val = qt90$pollination$per_group$`Medium.50%`[,1], lab = "Bees", col = cols_pol[1], type = "group"),
  # list(val = qt90$pollination$per_group$`Low.25%`[,1], lab = "L", col = cols_pol[1], type = "group"),
  list(type = "spacer", height = 0.5),
  
  #Moths
  # list(val = qt90$pollination$per_group$`High.75%`[,2], lab = "H", col = cols_pol[2], type = "group"),
  list(val = qt90$pollination$per_group$`Medium.50%`[,2], lab = "Moths", col = cols_pol[2], type = "group"),
  # list(val = qt90$pollination$per_group$`Low.25%`[,2], lab = "L", col = cols_pol[2], type = "group"),
  list(type = "spacer", height = 0.5),
  
  # Bats
  list(val = qt90$pollination$per_group$`High.75%`[,3], lab = "H", col = cols_pol[3], type = "group"),
  list(val = qt90$pollination$per_group$`Medium.50%`[,3], lab = "Bats", col = cols_pol[3], type = "group"),
  list(val = qt90$pollination$per_group$`Low.25%`[,3], lab = "L", col = cols_pol[3], type = "group"),
  list(type = "spacer", height = 0.5),
  
  #Multi
  # list(val = qt90$pollination$multi$`High.75%`, lab = "H", col = cols_pol[4], type ="multi"), 
  list(val = qt90$pollination$multi$`Medium.50%`, lab = "Multifunctionality", col = cols_pol[4], type ="multi"),
  # list(val = qt90$pollination$multi$`Low.25%`, lab = "L", col = cols_pol[4], type ="multi"),
  
  # Visual space
  list(type = "spacer", height = 2),
  
  # Seed dispersal
  
  #Bats
  list(val = qt90$dispersal$per_group$`High.75%`[,1], lab = "H", col = cols_sd[1], type = "group"),
  list(val = qt90$dispersal$per_group$`Medium.50%`[,1], lab = "M", col = cols_sd[1], type = "group"),
  list(val = qt90$dispersal$per_group$`Low.25%`[,1], lab = "L", col = cols_sd[1], type = "group"),
  list(type = "spacer", height = 0.5),
  
  # Birds
  list(val = qt90$dispersal$per_group$`High.75%`[,2], lab = "H", col = cols_sd[2], type = "group"),
  list(val = qt90$dispersal$per_group$`Medium.50%`[,2], lab = "M", col = cols_sd[2], type = "group"),
  list(val = qt90$dispersal$per_group$`Low.25%`[,2], lab = "L", col = cols_sd[2], type = "group"),
  list(type = "spacer", height = 0.5),
  
  # Non-flying mammals
  list(val = qt90$dispersal$per_group$`High.75%`[,3], lab = "H", col = cols_sd[3], type = "group"),
  list(val = qt90$dispersal$per_group$`Medium.50%`[,3], lab = "M", col = cols_sd[3], type = "group"),
  list(val = qt90$dispersal$per_group$`Low.25%`[,3], lab = "L", col = cols_sd[3], type = "group"),
  list(type = "spacer", height = 0.5),
  
  #Multi
  list(val = qt90$dispersal$multi$`High.75%`, lab = "H", col = cols_sd[4], type ="multi"),
  list(val = qt90$dispersal$multi$`Medium.50%`, lab = "M", col = cols_sd[4], type ="multi"),
  list(val = qt90$dispersal$multi$`Low.25%`, lab = "L", col = cols_sd[4], type ="multi")
)

# Definir limite do eixo X com base no maior valor de T90 observado nesses grupos
x_max <-max(as.numeric(qt90$dispersal$per_group$`Low.25%`)) * 1.1
# n_items <- length(summary_list) - 2

total_height <- sum(sapply(summary_list, function(x) if(!is.null(x$height)) x$height else 1))

plot(NULL, xlim = c(0, x_max), ylim = c(0.5, total_height + 0.5), 
     xlab = "Time to convergence (years)", ylab = "", 
     yaxt = "n", bty = "l")
mtext("C)", side = 3, line = 1.5, adj = 0, font = 2, cex = 1.2)

y_tracker <- total_height + 0.5

for(i in seq_along(summary_list)) {
  
  item <- summary_list[[i]]
  
  h <- if(!is.null(item$height)) item$height else 1
  
  y_pos <- y_tracker - (h / 2)
  y_tracker <- y_tracker - h
  
  if(is.null(item$type) || item$type == "spacer") next # Pular o espaço
  
  stats <- item$val
  lwd_thin <- ifelse(item$type == "multi", 2, 1)
  lwd_thick <- ifelse(item$type == "multi", 5, 3)
  
  # Barras de Erro
  segments(stats[1], y_pos, stats[5], y_pos, col = item$col, lwd = 2) # 90% CI
  segments(stats[2], y_pos, stats[4], y_pos, col = item$col, lwd = 6) # IQR
  points(stats[3], y_pos, pch = 21, bg = item$col, col = "white", cex = 2.8)
  
  # Nomes no Eixo Y
  font_type <- ifelse(item$type == "multi", 2, 1)
  # axis(2, at = y_pos, labels = item$lab, las = 1, cex.axis = 1, font = font_type, tick = FALSE)
  axis(2, at = y_pos, labels = item$lab, las = 1, cex.axis = 0.8, 
       tick = TRUE, tck = -0.02)
}

# Linha divisória entre os dois processos
# abline(h = 5, col = "grey80", lty = 2)

###################################################################################################################################
## Without connectivity
qt90 <- readRDS("output/t90.rds")
samples <- readRDS("output/model_posteriors.rds")
mcmc_mat <- as.matrix(samples)

weights_df <- read.csv("data/processed/weights_df.csv")

data <- read.csv("data/processed/model_df.csv")
data$type <- factor(ifelse(1:nrow(data) %in% grep("OG", data$Plot_ID), "old", "rec"),
                    levels = c("old", "rec"))


data <- data %>% mutate(ConIndex = as.numeric(scale(ConIndex)))

rec_conn <- data$ConIndex[data$type == "rec"]
conn_values <- c(
  Low    = quantile(rec_conn, 0.25, na.rm = TRUE),
  Medium = quantile(rec_conn, 0.50, na.rm = TRUE),
  High   = quantile(rec_conn, 0.75, na.rm = TRUE)
)

dataSub <- subset(data, select = c(type, RegTime,  
                                   FDBees, FDMoths,  FDBat_pol, 
                                   FDBats, FDBirds, FDNf, 
                                   FDSeeds, AbSeeds, RichSeeds, ShnSeeds,
                                   FDSdlng, AbSdlng, RichSdlng, ShnSdlng)) 
dataSub <- dataSub %>%
  mutate(across(
    .cols = -c(1:2),     
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
  pivot_longer(cols = -c(1:2), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = colnames(dataSub)[-c(1:2)]), variable = as.integer(variable)) %>%
  filter(!is.na(value))

### Recovery trajectories

var_names <- colnames(dataSub)[-c(1,2)]

#### Multifunctionality

cols_pol <- c("#cbc9e2", "#9e9ac8","#c994c7", "#3f007d")
cols_sd <- c("#3182bd", "#7bccc4", "#bdd7e7", "#08306b")

# cols_pol <- c("#3f007d", "#807dba", "#c994c7")
# cols_sd <- c("#08306b", "#3182bd", "#9ecae1")

# Calculate Weighted Data Points:

w_pol <- c(FDBees=weights_df[1,"prctg"]/sum(weights_df[1:3,"prctg"]), FDMoths=weights_df[2,"prctg"]/sum(weights_df[1:3,"prctg"]), FDBat_pol=weights_df[3,"prctg"]/sum(weights_df[1:3,"prctg"])) 
w_sd <- c(FDBats=weights_df[4,"prctg"]/sum(weights_df[4:6,"prctg"]), FDBirds=weights_df[5,"prctg"]/sum(weights_df[4:6,"prctg"]), FDNf=weights_df[6,"prctg"]/sum(weights_df[4:6,"prctg"]))

dataSub$multi_pol <- apply(dataSub[, names(w_pol)], 1, function(x) {
  sum(x * w_pol, na.rm = TRUE)
})
dataSub$multi_sd <-apply(dataSub[, names(w_sd)], 1, function(x) {
  sum(x * w_sd, na.rm = TRUE)
})

compute_multi_trajectory <- function(draw, idxs, weights, conn, t) {
  
  At <- numeric(length(idxs))
  
  for(i in seq_along(idxs)) {
    j <- idxs[i]
    
    t0   <- draw[paste0("theta_0[", j, "]")]
    tinf <- draw[paste0("theta_inf[", j, "]")]
    lambda <- draw[paste0("lambda[", j, "]")]
    
    At[i] <- t0 + (tinf - t0) * (1 - exp(-lambda * t))
  }
  
  sum(At * weights)
}

plot_trajectory_multi <- function(var_indices, weights,
                                  conn, title, colors,
                                  data_subset) {
  
  plot(NULL, xlim = c(0, 40), ylim = c(0, 1.2), 
       xlab = "Time (years)", ylab = "Functional Diversity",
       bty = "l")
  mtext(title, side = 3, line = 1.5, adj = 0, font = 2, cex = 1.2)
  
  # ----- 1) Individual groups (as you already had) -----
  
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
            col = adjustcolor(colors[i], alpha.f = 0.1),
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
    
    points(data_subset$RegTime[data_subset$variable == idx], 
           data_subset$value[data_subset$variable == idx], 
           pch = 21,
           bg = adjustcolor(colors[i], alpha.f = 1),
           col = "white", cex = 1.2)
  }
  
  # ----- 2) MULTIFUNCTIONALITY -----
  
  # t_seq <- seq(0, 40, length = 200)
  
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
  
  # Median line multifunctionality
  lines(t_seq, med, col = colors[4], lwd = 4)
  
  # legend("bottomright",
  #        legend = c(names(weights), "Multifunctionality"),
  #        col = c(colors, "black"),
  #        lwd = c(rep(3, length(colors)), 4),
  #        bty = "n")
}



pdf("Figure_trends.pdf", width = 12, height = 5)

layout_mat <- matrix(c(1, 2, 3, 4, 5), nrow = 1)
layout(layout_mat, widths = c(3, 0.8, 3, 0.8, 4.4))

par(mar = c(4, 4, 3, 0.8), oma = c(1, 1, 1, 1), mgp = c(2, 0.7, 0), tcl = -0.3)

plot_old_growth_multi <- function(var_indices, weights, colors) {
  
  n <- length(var_indices)
  
  plot(NULL,
       xlim = c(0.5, n + 1.5),
       ylim = c(0, 1.2), 
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
  
  # Labels
  # labs <- c(names(weights), "Multifunctionality")
  
  axis(1,
       at = 1:(n+1),
       labels = FALSE,
       lwd.ticks = 0,
       las = 2, cex.axis = 0.8)
  
  mtext("Old-growth", side = 1, line = 2.1, cex = 0.7)
}
# --- PAINEL A: POLLINATION ---
par(mar = c(4, 4, 3, 0.8))
plot_trajectory_multi(1:3, w_pol, conn_values["Medium.50%"], "A) Pollination", cols_pol, long[long$type == "rec", ])
par(mar = c(4, 0.8, 3, 1))
plot_old_growth_multi(1:3, w_pol, cols_pol)

# --- PAINEL B: SEED DISPERSAL ---
par(mar = c(4, 4, 3, 0.8))
plot_trajectory_multi(4:6, w_sd, conn_values["Medium.50%"], "B) Seed dispersal", cols_sd, long[long$type == "rec", ])
par(mar = c(4, 0.8, 3, 1))
plot_old_growth_multi(4:6, w_sd, cols_sd)

# --- PAINEL C: RECOVERY SUMMARY (T90) ---
par(mar = c(4, 10, 3, 1))

# summary_list <- list(
#   
#   # Pollination
#   
#   #High
#   list(val = qt90$pollination$multi$`High.75%`, lab = "", col = cols_pol[4], type ="multi"),
#   list(val = qt90$pollination$per_group$`High.75%`[,1], lab = "", col = cols_pol[1], type = "group"),
#   list(val = qt90$pollination$per_group$`High.75%`[,2], lab = "", col = cols_pol[2], type = "group"),
#   list(val = qt90$pollination$per_group$`High.75%`[,3], lab = "", col = cols_pol[3], type = "group"),
#   list(type = "spacer", height = 0.5),
#   #Medium
#   list(val = qt90$pollination$multi$`Medium.50%`, lab = "", col = cols_pol[4], type ="multi"),
#   list(val = qt90$pollination$per_group$`Medium.50%`[,1], lab = "", col = cols_pol[1], type = "group"),
#   list(val = qt90$pollination$per_group$`Medium.50%`[,2], lab = "", col = cols_pol[2], type = "group"),
#   list(val = qt90$pollination$per_group$`Medium.50%`[,3], lab = "", col = cols_pol[3], type = "group"),
#   list(type = "spacer", height = 0.5),
#   #Low
#   list(val = qt90$pollination$multi$`Low.25%`, lab = "", col = cols_pol[4], type ="multi"),
#   list(val = qt90$pollination$per_group$`Low.25%`[,1], lab = "", col = cols_pol[1], type = "group"),
#   list(val = qt90$pollination$per_group$`Low.25%`[,2], lab = "", col = cols_pol[2], type = "group"),
#   list(val = qt90$pollination$per_group$`Low.25%`[,3], lab = "", col = cols_pol[3], type = "group"),
#   
#   # Visual space
#   list(type = "spacer", height = 2),
#   
#   # Seed dispersal
#   
#   #High
#   list(val = qt90$dispersal$multi$`High.75%`, lab = "", col = cols_sd[4], type ="multi"),
#   list(val = qt90$dispersal$per_group$`High.75%`[,1], lab = "", col = cols_sd[1], type = "group"),
#   list(val = qt90$dispersal$per_group$`High.75%`[,2], lab = "", col = cols_sd[2], type = "group"),
#   list(val = qt90$dispersal$per_group$`High.75%`[,3], lab = "", col = cols_sd[3], type = "group"),
#   list(type = "spacer", height = 0.5),
#   #Medium
#   list(val = qt90$dispersal$multi$`Medium.50%`, lab = "", col = cols_sd[4], type ="multi"),
#   list(val = qt90$dispersal$per_group$`Medium.50%`[,1], lab = "", col = cols_sd[1], type = "group"),
#   list(val = qt90$dispersal$per_group$`Medium.50%`[,2], lab = "", col = cols_sd[2], type = "group"),
#   list(val = qt90$dispersal$per_group$`Medium.50%`[,3], lab = "", col = cols_sd[3], type = "group"),
#   list(type = "spacer", height = 0.5),
#   #Low
#   list(val = qt90$dispersal$multi$`Low.25%`, lab = "", col = cols_sd[4], type ="group"),
#   list(val = qt90$dispersal$per_group$`Low.25%`[,1], lab = "", col = cols_sd[1], type = "group"),
#   list(val = qt90$dispersal$per_group$`Low.25%`[,2], lab = "", col = cols_sd[2], type = "group"),
#   list(val = qt90$dispersal$per_group$`Low.25%`[,3], lab = "", col = cols_sd[3], type = "group")
# )

summary_list <- list(
  
  # Pollination
  
  #Multi
  list(val = qt90$pollination$multi$`High.75%`, lab = "H", col = cols_pol[4], type ="multi"), 
  list(val = qt90$pollination$multi$`Medium.50%`, lab = "M", col = cols_pol[4], type ="multi"),
  list(val = qt90$pollination$multi$`Low.25%`, lab = "L", col = cols_pol[4], type ="multi"),
  list(type = "spacer", height = 0.5),
  
  #Bees
  list(val = qt90$pollination$per_group$`High.75%`[,1], lab = "H", col = cols_pol[1], type = "group"),
  list(val = qt90$pollination$per_group$`Medium.50%`[,1], lab = "M", col = cols_pol[1], type = "group"),
  list(val = qt90$pollination$per_group$`Low.25%`[,1], lab = "L", col = cols_pol[1], type = "group"),
  list(type = "spacer", height = 0.5),
  
  #Moths
  list(val = qt90$pollination$per_group$`High.75%`[,2], lab = "H", col = cols_pol[2], type = "group"),
  list(val = qt90$pollination$per_group$`Medium.50%`[,2], lab = "M", col = cols_pol[2], type = "group"),
  list(val = qt90$pollination$per_group$`Low.25%`[,2], lab = "L", col = cols_pol[2], type = "group"),
  list(type = "spacer", height = 0.5),
  
  # Bats
  list(val = qt90$pollination$per_group$`High.75%`[,3], lab = "H", col = cols_pol[3], type = "group"),
  list(val = qt90$pollination$per_group$`Medium.50%`[,3], lab = "M", col = cols_pol[3], type = "group"),
  list(val = qt90$pollination$per_group$`Low.25%`[,3], lab = "L", col = cols_pol[3], type = "group"),
  
  # Visual space
  list(type = "spacer", height = 2),
  
  # Seed dispersal
  
  #Multi
  list(val = qt90$dispersal$multi$`High.75%`, lab = "H", col = cols_sd[4], type ="multi"),
  list(val = qt90$dispersal$multi$`Medium.50%`, lab = "M", col = cols_sd[4], type ="multi"),
  list(val = qt90$dispersal$multi$`Low.25%`, lab = "L", col = cols_sd[4], type ="multi"),
  list(type = "spacer", height = 0.5),
  
  #Bats
  list(val = qt90$dispersal$per_group$`High.75%`[,1], lab = "H", col = cols_sd[1], type = "group"),
  list(val = qt90$dispersal$per_group$`Medium.50%`[,1], lab = "M", col = cols_sd[1], type = "group"),
  list(val = qt90$dispersal$per_group$`Low.25%`[,1], lab = "L", col = cols_sd[1], type = "group"),
  list(type = "spacer", height = 0.5),
  
  # Birds
  list(val = qt90$dispersal$per_group$`High.75%`[,2], lab = "H", col = cols_sd[2], type = "group"),
  list(val = qt90$dispersal$per_group$`Medium.50%`[,2], lab = "M", col = cols_sd[2], type = "group"),
  list(val = qt90$dispersal$per_group$`Low.25%`[,2], lab = "L", col = cols_sd[2], type = "group"),
  list(type = "spacer", height = 0.5),
  
  # Non-flying mammals
  list(val = qt90$dispersal$per_group$`High.75%`[,3], lab = "H", col = cols_sd[3], type = "group"),
  list(val = qt90$dispersal$per_group$`Medium.50%`[,3], lab = "M", col = cols_sd[3], type = "group"),
  list(val = qt90$dispersal$per_group$`Low.25%`[,3], lab = "L", col = cols_sd[3], type = "group")
)

# Definir limite do eixo X com base no maior valor de T90 observado nesses grupos
x_max <-max(as.numeric(qt90$dispersal$per_group$`Low.25%`)) * 1.1
# n_items <- length(summary_list) - 2

total_height <- sum(sapply(summary_list, function(x) if(!is.null(x$height)) x$height else 1))

plot(NULL, xlim = c(0, x_max), ylim = c(0.5, total_height + 0.5), 
     xlab = "Time to convergence (years)", ylab = "", 
     yaxt = "n", bty = "l")
mtext("C)", side = 3, line = 1.5, adj = 0, font = 2, cex = 1.2)

y_tracker <- total_height + 0.5

for(i in seq_along(summary_list)) {
  
  item <- summary_list[[i]]
  
  h <- if(!is.null(item$height)) item$height else 1
  
  y_pos <- y_tracker - (h / 2)
  y_tracker <- y_tracker - h
  
  if(is.null(item$type) || item$type == "spacer") next # Pular o espaço
  
  stats <- item$val
  lwd_thin <- ifelse(item$type == "multi", 2, 1)
  lwd_thick <- ifelse(item$type == "multi", 5, 3)
  
  # Barras de Erro
  segments(stats[1], y_pos, stats[5], y_pos, col = item$col, lwd = 2) # 90% CI
  segments(stats[2], y_pos, stats[4], y_pos, col = item$col, lwd = 6) # IQR
  points(stats[3], y_pos, pch = 21, bg = item$col, col = "white", cex = 2.8)
  
  # Nomes no Eixo Y
  font_type <- ifelse(item$type == "multi", 2, 1)
  # axis(2, at = y_pos, labels = item$lab, las = 1, cex.axis = 1, font = font_type, tick = FALSE)
  axis(2, at = y_pos, labels = item$lab, las = 1, cex.axis = 0.8, 
       tick = TRUE, tck = -0.02)
}

# Linha divisória entre os dois processos
# abline(h = 5, col = "grey80", lty = 2)

dev.off()
