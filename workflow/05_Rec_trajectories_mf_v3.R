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
    a    <- draw[paste0("alpha_con[", j, "]")]
    b    <- draw[paste0("beta_con[", j, "]")]
    
    lambda <- exp(a + b * conn)
    
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
plot_trajectory_multi(1:3, w_pol, conn_values["Medium.50%"], "A)", cols_pol, long[long$type == "rec", ])
par(mar = c(4, 0.8, 3, 1))
plot_old_growth_multi(1:3, w_pol, cols_pol)

# --- PAINEL B: SEED DISPERSAL ---
par(mar = c(4, 4, 3, 0.8))
plot_trajectory_multi(4:6, w_sd, conn_values["Medium.50%"], "B)", cols_sd, long[long$type == "rec", ])
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
     xlab = "Time to old-growth (years)", ylab = "", 
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
