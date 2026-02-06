#### Data

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

#### Figure 1: Multifunctionality

# Colors
cols_multi <- c("#3f007d", "#1F78B4") 

# Calculate Weighted Data Points:

w_pol <- c(FDBees=weights_df[1,"prctg"]/sum(weights_df[1:3,"prctg"]), FDMoths=weights_df[2,"prctg"]/sum(weights_df[1:3,"prctg"]), FDBat_pol=weights_df[3,"prctg"]/sum(weights_df[1:3,"prctg"])) 
w_sd <- c(FDBats=weights_df[4,"prctg"]/sum(weights_df[4:6,"prctg"]), FDBirds=weights_df[5,"prctg"]/sum(weights_df[4:6,"prctg"]), FDNf=weights_df[6,"prctg"]/sum(weights_df[4:6,"prctg"]))

dataSub$multi_pol <- apply(dataSub[, names(w_pol)], 1, function(x) {
  sum(x * w_pol, na.rm = TRUE)
})
dataSub$multi_sd <-apply(dataSub[, names(w_sd)], 1, function(x) {
  sum(x * w_sd, na.rm = TRUE)
})

get_multi_recovery_traj <- function(samples, groups, conn, weights, t_seq) {
  
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
  
  lambda <- exp(aMat + bMat * conn)
  
  traj <- matrix(NA, nrow(theta_0), length(t_seq))
  
  for(j in seq_along(t_seq)) {
    tt <- t_seq[j]
    
    At <- theta_0 + (theta_inf - theta_0) * (1 - exp(-lambda * tt))
    dev_rel <- abs(At - theta_inf) / theta_inf
    Mt <- as.numeric(dev_rel %*% w)
    
    traj[, j] <- 1 - Mt   # RECOVERY %
  }
  
  traj
}

layout(matrix(c(1, 2), nrow = 1), widths = c(1.2, 1))
par(oma = c(2, 1, 1, 1), mar = c(4, 4, 3, 1), las = 1, tcl = -0.3)

# --- PANEL A: Combined Trajectories ---
plot(NULL, xlim = c(0, 40), ylim = c(0, 1), bty = "l",
     xlab = "Time (years)", ylab = "Relative recovery", 
     main = "A) Recovery Trajectories")

draw_recovery_traj <- function(traj, col, t_seq, show_threshold = TRUE) {
  
  polygon(
    c(t_seq, rev(t_seq)),
    c(apply(traj, 2, quantile, 0.025),
      rev(apply(traj, 2, quantile, 0.975))),
    col = adjustcolor(col, 0.15),
    border = NA
  )
  
  lines(t_seq, apply(traj, 2, median), col = col, lwd = 3)
  
  if(show_threshold) {
    abline(h = 0.9, lty = 2, col = "grey40")
  }
}

t_seq <- seq(0, 40, length.out = 100)

# Pollination
traj_pol <- get_multi_recovery_traj(
  samples,
  groups = 1:3,
  conn   = conn_values["Medium.50%"],
  weights = w_pol,
  t_seq = t_seq
)
draw_recovery_traj(traj_pol, cols_multi[1], t_seq)

# Seed dispersal
traj_sd <- get_multi_recovery_traj(
  samples,
  groups = 4:6,
  conn   = conn_values["Medium.50%"],
  weights = w_sd,
  t_seq = t_seq
)
draw_recovery_traj(traj_sd, cols_multi[2], t_seq)

legend("bottomright", legend = c("Pollination", "Seed Dispersal"), 
       col = cols_multi, lwd = 3, bty = "n", cex = 0.9)

# --- PANEL B: Recovery Time ---
par(mar = c(4, 5, 3, 2)) 

plot(NULL, xlim = c(0, 58.8*1.1), ylim = c(0.5, 6.5), yaxt = "n", 
     xlab = "Years to recovery", ylab = "", main = "B) Time to Recovery", bty = "l")

procs_to_plot <- c("pollination", "dispersal")
levels_to_plot <- c("High.75%", "Medium.50%", "Low.25%")
y_val <- 6

for(proc in procs_to_plot) {
  # # Add section headers for Connectivity
  # header_pos <- if(proc == "pollination") 6.4 else 3.4
  # mtext(ifelse(proc == "pollination", "Pollination", "Seed Dispersal"), 
  #       side = 2, at = header_pos, line = 3, las = 0, font = 2, cex = 0.8, adj = 1)
  # 
  for(lvl in levels_to_plot) {
    curr_col <- if(proc == "pollination") cols_multi[1] else cols_multi[2]
    stats <- qt90[[proc]]$multi[[lvl]]
    
    segments(stats[1], y_val, stats[5], y_val, col = curr_col, lwd = 1.5)
    segments(stats[2], y_val, stats[4], y_val, col = curr_col, lwd = 5.5)
    points(stats[3], y_val, pch = 21, bg = "white", col = curr_col, cex = 1.5, lwd = 2)
    
    # Labeling: Only High, Medium, Low
    clean_lvl <- gsub("\\.[0-9]+%", "", lvl)
    axis(2, at = y_val, labels = clean_lvl, cex.axis = 0.8)
    
    y_val <- y_val - 1
  }
}

#### Figure 2: Groups

cols_pol <- c("#3f007d", "#807dba", "#c994c7")
cols_sd <- c("#08306b", "#3182bd", "#9ecae1")

# pdf("recovery_interactions_only.pdf", width = 12, height = 5)

layout_mat <- matrix(c(1, 2, 3, 4, 5), nrow = 1)
layout(layout_mat, widths = c(3, 0.8, 3, 0.8, 4.4))

par(mar = c(4, 4, 3, 0.8), oma = c(1, 1, 1, 1), mgp = c(2, 0.7, 0), tcl = -0.3)

plot_trajectory <- function(var_indices, conn, title, colors, data_subset) {
  plot(NULL, xlim = c(0, 40), ylim = c(0, 1.2), 
       xlab = "Time (years)", ylab = "Functional Diversity", main = title, bty = "l")
  
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
            add = TRUE, col = adjustcolor(colors[i], alpha.f = 0.1), lwd = 0.5)
    }
    
    # Mediana
    t0_m <- median(mcmc_mat[, paste0("theta_0[", idx, "]")])
    tinf_m <- median(mcmc_mat[, paste0("theta_inf[", idx, "]")])
    a_m <- median(mcmc_mat[, paste0("alpha_con[", idx, "]")])
    b_m <- median(mcmc_mat[, paste0("beta_con[", idx, "]")])
    
    lamb_m <- exp(a_m + b_m * conn)
    
    curve(t0_m + (tinf_m - t0_m) * (1 - exp(-lamb_m * x)), 
          add = TRUE, col = colors[i], lwd = 3)
    
    points(data_subset$RegTime[data_subset$variable == idx], 
           data_subset$value[data_subset$variable == idx], 
           pch = 21, bg = adjustcolor(colors[i], alpha.f = 1), col = "white", cex = 1.2)
  }
}

plot_old_growth <- function(var_indices, colors) {
  plot(NULL, xlim = c(0.5, length(var_indices) + 0.5), ylim = c(0, 1.2), 
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  # abline(h = seq(0, 1, 0.2), col = "grey90", lty = 2)
  axis(1, at = c(0.5, length(var_indices) + 0.5), labels = FALSE, lwd.ticks = 0, lwd = 1)
  
  for(i in seq_along(var_indices)) {
    idx <- var_indices[i]
    vals <- mcmc_mat[, paste0("theta_inf[", idx, "]")]
    stats <- quantile(vals, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    
    # Desenhar barras de erro
    segments(i, stats[1], i, stats[5], col = colors[i], lwd = 1)
    segments(i, stats[2], i, stats[4], col = colors[i], lwd = 4)
    points(i, stats[3], pch = 21, bg = colors[i], col = colors[i], cex = 1.5, lwd = 2)
  }
  mtext("Old-growth", side = 1, line = 0.5, cex = 0.7)
}
# --- PAINEL A: POLLINATION ---
par(mar = c(4, 4, 3, 0.8))
plot_trajectory(1:3, conn_values["Medium.50%"], "A) Pollination", cols_pol, long[long$type == "rec", ])
par(mar = c(4, 0.8, 3, 1))
plot_old_growth(1:3, cols_pol)

# --- PAINEL B: SEED DISPERSAL ---
par(mar = c(4, 4, 3, 0.8))
plot_trajectory(4:6, conn_values["Medium.50%"], "B) Seed Dispersal", cols_sd, long[long$type == "rec", ])
par(mar = c(4, 0.8, 3, 1))
plot_old_growth(4:6, cols_sd)

# --- PAINEL C: RECOVERY SUMMARY (T90) ---
par(mar = c(4, 10, 3, 1))
n_items <- length(summary_list)

summary_list <- list(
  
  # Pollination
  list(val = qt90$pollination$per_group$`Medium.50%`[,1], lab = "Bees", col = cols_pol[1], type = "group"),
  list(val = qt90$pollination$per_group$`Medium.50%`[,2], lab = "Moths", col = cols_pol[2], type = "group"),
  list(val = qt90$pollination$per_group$`Medium.50%`[,3], lab = "Bats", col = cols_pol[3], type = "group"),
  
  # Visual space
  list(type = "spacer"),
  
  # Seed dispersal
  list(val = qt90$dispersal$per_group$`Medium.50%`[,1], lab = "Bats", col = cols_sd[1], type = "group"),
  list(val = qt90$dispersal$per_group$`Medium.50%`[,2], lab = "Birds", col = cols_sd[2], type = "group"),
  list(val = qt90$dispersal$per_group$`Medium.50%`[,3], lab = "Non-flying mammals", col = cols_sd[3], type = "group")
)

# Definir limite do eixo X com base no maior valor de T90 observado nesses grupos
x_max <-max(as.numeric(qt90$dispersal$per_group$`Medium.50%`)) * 1.1
n_items <- 7

plot(NULL, xlim = c(0, x_max), ylim = c(0.5, n_items + 0.5), 
     xlab = "Years to recovery", ylab = "", 
     yaxt = "n", main = "C) Time to Recovery", bty = "l")

for(i in 1:n_items) {
  item <- summary_list[[i]]
  y_pos <- n_items - i + 1
  
  if(is.null(item$type) || item$type == "spacer") next # Pular o espaço
  
  stats <- item$val
  # lwd_thin <- ifelse(item$type == "multi", 2, 1)
  # lwd_thick <- ifelse(item$type == "multi", 5, 3)
  
  # Barras de Erro
  segments(stats[1], y_pos, stats[5], y_pos, col = item$col, lwd = 2) # 90% CI
  segments(stats[2], y_pos, stats[4], y_pos, col = item$col, lwd = 6) # IQR
  points(stats[3], y_pos, pch = 21, bg = item$col, col = "white", cex = 2.8)
  
  # Nomes no Eixo Y
  font_type <- ifelse(item$type == "multi", 2, 1)
  axis(2, at = y_pos, labels = item$lab, las = 1, cex.axis = 1, font = font_type, tick = FALSE)
}

# Linha divisória entre os dois processos
# abline(h = 5, col = "grey80", lty = 2)

# dev.off()
