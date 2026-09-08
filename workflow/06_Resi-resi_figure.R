#### Resistance and Resilience Figure ####
#
# Resistance: abs(1 - (theta_0 - theta_inf) / theta_inf) per posterior draw.
#             Values close to 1 = high resistance (close to old-growth at t=0).
#             Values far from 1 = low resistance (strongly deviated from old-growth).
#
# Resilience: lambda at median connectivity = exp(alpha_con + beta_con * conn_median).
#             Higher lambda = faster recovery.
#
# No new model runs needed — all quantities derivable from saved posteriors.

#### Data ####

samples  <- readRDS("output/model_posteriors_con.rds")
mcmc_mat <- as.matrix(samples)

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

#### Colors ####

cols_pol <- c("#c994c7", "#9e9ac8", "#cbc9e2", "#3f007d")  # Bats, Bees, Moths, multi
cols_sd  <- c("#3182bd", "#7bccc4", "#bdd7e7", "#08306b")  # Bats, Birds, NF mammals, multi

# group indices: 1=Bees, 2=Moths, 3=Bat_pol, 4=Bat_sd, 5=Birds, 6=NF

#### Posterior extraction helpers ####

get_resistance <- function(j) {
  theta_0   <- mcmc_mat[, paste0("theta_0[",   j, "]")]
  theta_inf <- mcmc_mat[, paste0("theta_inf[", j, "]")]
  1 / (1 + abs(theta_0 - theta_inf))
}

get_resilience <- function(j, conn) {
  alpha <- mcmc_mat[, paste0("alpha_con[", j, "]")]
  beta  <- mcmc_mat[, paste0("beta_con[",  j, "]")]
  exp(alpha + beta * conn)
}

ci_summary <- function(x) {
  quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
}

#### Build summary lists ####

conn_med <- conn_values["Medium.50%"]

# --- Pollination (Bats=3, Bees=1, Moths=2) ---

resistance_pol <- list(
  list(val = ci_summary(get_resistance(3)), lab = "Bats",  col = cols_pol[1], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = ci_summary(get_resistance(1)), lab = "Bees",  col = cols_pol[2], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = ci_summary(get_resistance(2)), lab = "Moths", col = cols_pol[3], type = "group")
)

resilience_pol <- list(
  list(val = ci_summary(get_resilience(3, conn_med)), lab = "Bats",  col = cols_pol[1], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = ci_summary(get_resilience(1, conn_med)), lab = "Bees",  col = cols_pol[2], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = ci_summary(get_resilience(2, conn_med)), lab = "Moths", col = cols_pol[3], type = "group")
)

# --- Seed dispersal (Bats=4, Birds=5, NF=6) ---

resistance_sd <- list(
  list(val = ci_summary(get_resistance(4)), lab = "Bats",               col = cols_sd[1], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = ci_summary(get_resistance(5)), lab = "Birds",              col = cols_sd[2], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = ci_summary(get_resistance(6)), lab = "Non-flying\nmammals", col = cols_sd[3], type = "group")
)

resilience_sd <- list(
  list(val = ci_summary(get_resilience(4, conn_med)), lab = "Bats",               col = cols_sd[1], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = ci_summary(get_resilience(5, conn_med)), lab = "Birds",              col = cols_sd[2], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = ci_summary(get_resilience(6, conn_med)), lab = "Non-flying\nmammals", col = cols_sd[3], type = "group")
)

#### Plotting function (matches draw_rectime_plot style) ####

draw_stability_plot <- function(current_list, label, x_lim, xlab) {
  
  total_height <- sum(sapply(current_list, function(x) if (!is.null(x$height)) x$height else 1))
  
  plot(NULL, xlim = x_lim, ylim = c(0.5, total_height + 0.5),
       xlab = xlab, ylab = "",
       yaxt = "n", bty = "l")
  
  if (!is.null(label)) mtext(label, side = 3, line = 1.5, adj = 0, font = 2, cex = 1.2)
  
  y_tracker <- total_height + 0.5
  
  for (i in seq_along(current_list)) {
    item  <- current_list[[i]]
    h     <- if (!is.null(item$height)) item$height else 1
    y_pos <- y_tracker - (h / 2)
    y_tracker <- y_tracker - h
    
    if (is.null(item$type) || item$type == "spacer") next
    
    stats <- item$val
    
    segments(stats[1], y_pos, stats[5], y_pos, col = item$col, lwd = 2)  # 90% CI
    segments(stats[2], y_pos, stats[4], y_pos, col = item$col, lwd = 6)  # IQR
    points(stats[3], y_pos, pch = 21, bg = item$col, col = "white", cex = 2.5)
    
    axis(2, at = y_pos, labels = item$lab, las = 1, cex.axis = 0.9,
         font = 1, tick = TRUE, tck = -0.02)
  }
}

#### X limits (per panel, from data) ####

x_max <- function(item_list) {
  max(unlist(lapply(item_list, function(x) if (!is.null(x$val)) x$val else NULL)),
      na.rm = TRUE) * 1.1
}

res_pol_x_lim <- c(0.5, 1)
res_sd_x_lim  <- c(0.5, 1)
resil_pol_x_max <- x_max(resilience_pol)
resil_sd_x_max  <- x_max(resilience_sd)

#### Figure ####

svg("output/Figures/Figure_resistance_resilience.svg", width = 10, height = 4)

layout(matrix(c(1, 2, 3, 4), nrow = 1), widths = c(3, 2, 3, 2))
par(oma = c(1, 1, 1, 1), mgp = c(2, 0.7, 0), tcl = -0.3)

# --- Panel A: Pollination resistance ---
par(mar = c(4, 6, 3, 1))
draw_stability_plot(resistance_pol, "A) Pollination", c(0.5, 1), "Resistance")

# --- Panel B: Pollination resilience ---
par(mar = c(4, 4, 3, 1))
draw_stability_plot(resilience_pol, NULL, c(0, resil_pol_x_max), "Resilience")

# --- Panel C: Seed dispersal resistance ---
par(mar = c(4, 6, 3, 1))
draw_stability_plot(resistance_sd, "B) Seed dispersal", c(0.5, 1), "Resistance")

# --- Panel D: Seed dispersal resilience ---
par(mar = c(4, 4, 3, 1))
draw_stability_plot(resilience_sd, NULL, c(0, resil_sd_x_max), "Resilience")

dev.off()

######## Sanity check

metadata <- readRDS("output/model_metadata_con.rds")
metadata$var_names
colnames(dataSub)[-c(1:2)]

theta_0_birds   <- median(mcmc_mat[, "theta_0[5]"])
theta_inf_birds <- median(mcmc_mat[, "theta_inf[5]"])

theta_0_birds
theta_inf_birds

# resistance as implemented
abs(1 - (theta_0_birds - theta_inf_birds) / theta_inf_birds)
1/(1 + abs(theta_0_birds - theta_inf_birds))

data %>% 
  filter(type == "rec") %>% 
  summarise(mean_birds = mean(FDBirds, na.rm = TRUE))

data %>% 
  filter(type == "old") %>% 
  summarise(mean_birds = mean(FDBirds, na.rm = TRUE))

grep("theta_inf", colnames(mcmc_mat), value = TRUE)

median(mcmc_mat[, "theta_inf[5]"])
median(mcmc_mat[, "theta_0[5]"])

median(mcmc_mat[, "theta_inf[4]"])
median(mcmc_mat[, "theta_0[4]"])

theta_0_bats   <- median(mcmc_mat[, "theta_0[4]"])
theta_inf_bats <- median(mcmc_mat[, "theta_inf[4]"])

abs(1 - (theta_0_bats - theta_inf_bats) / theta_inf_bats)
1/(1 + abs(theta_0_bats - theta_inf_bats))

theta_0_nf   <- median(mcmc_mat[, "theta_0[6]"])
theta_inf_nf <- median(mcmc_mat[, "theta_inf[6]"])

abs(1 - (theta_0_nf - theta_inf_nf) / theta_inf_nf)
1/(1 + abs(theta_0_nf - theta_inf_nf))

theta_0_bats   <- median(mcmc_mat[, "theta_0[4]"])
theta_inf_bats <- median(mcmc_mat[, "theta_inf[4]"])

abs(1 - (theta_0_bats - theta_inf_bats) / theta_inf_bats)
1/(1 + abs(theta_0_bats - theta_inf_bats))

theta_0_batspol   <- median(mcmc_mat[, "theta_0[1]"])
theta_inf_batspol <- median(mcmc_mat[, "theta_inf[1]"])

abs(1 - (theta_0_batspol - theta_inf_batspol) / theta_inf_batspol)
1/(1 + abs(theta_0_batspol - theta_inf_batspol))

1/(1 + abs(median(mcmc_mat[, "theta_0[2]"]) - median(mcmc_mat[, "theta_inf[2]"])))
1/(1 + abs(median(mcmc_mat[, "theta_0[3]"]) - median(mcmc_mat[, "theta_inf[3]"])))


#### Per-group net change in contribution during recovery ####
#
# Buffering is visible as offsetting net changes across groups:
# each group's weighted contribution changes from active land use (theta_0)
# to old-growth (theta_inf) by  w_k * (theta_inf_k - theta_0_k).
# Positive = group's contribution increases during recovery.
# Negative = group's contribution decreases during recovery.
# The multi-group net change is the sum across groups; near-zero = buffered.
#
# Computed per posterior draw (difference of two parameters, so far better
# constrained than variance-ratio metrics). No model rerun needed.

#### Data ####

samples  <- readRDS("output/model_posteriors_con.rds")
mcmc_mat <- as.matrix(samples)

weights_df <- read.csv("data/processed/weights_df.csv")

#### Colors ####

cols_pol <- c("#c994c7", "#9e9ac8", "#cbc9e2", "#3f007d")  # Bats, Bees, Moths, multi
cols_sd  <- c("#3182bd", "#7bccc4", "#bdd7e7", "#08306b")  # Bats, Birds, NF mammals, multi

#### Weights (same as trajectory figure) ####

w_pol <- c(FDBees  = weights_df[1,"prctg"]/sum(weights_df[1:3,"prctg"]),
           FDMoths = weights_df[2,"prctg"]/sum(weights_df[1:3,"prctg"]),
           FDBat_pol = weights_df[3,"prctg"]/sum(weights_df[1:3,"prctg"]))
w_sd  <- c(FDBats  = weights_df[4,"prctg"]/sum(weights_df[4:6,"prctg"]),
           FDBirds = weights_df[5,"prctg"]/sum(weights_df[4:6,"prctg"]),
           FDNf    = weights_df[6,"prctg"]/sum(weights_df[4:6,"prctg"]))

#### Posterior of net change ####

# Per-group weighted net change: w_k * (theta_inf_k - theta_0_k), per draw.
# get_net_change <- function(j, w) {
#   theta_0   <- mcmc_mat[, paste0("theta_0[",   j, "]")]
#   theta_inf <- mcmc_mat[, paste0("theta_inf[", j, "]")]
#   w * (theta_inf - theta_0)
# }
get_net_change <- function(j) {
  theta_0   <- mcmc_mat[, paste0("theta_0[",   j, "]")]
  theta_inf <- mcmc_mat[, paste0("theta_inf[", j, "]")]
  theta_inf - theta_0
}

# Multi-group net change: sum of weighted group net changes, per draw.
# get_multi_net_change <- function(idxs, weights) {
#   out <- 0
#   for (i in seq_along(idxs)) {
#     out <- out + get_net_change(idxs[i], weights[i])
#   }
#   out
# }
get_multi_net_change <- function(idxs, weights) {
  out <- numeric(nrow(mcmc_mat))
  for (i in seq_along(idxs)) {
    out <- out + weights[i] * get_net_change(idxs[i])
  }
  out
}

ci_summary <- function(x) quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)

#### Build summary lists (plotting order matches trajectory figure) ####

# Pollination: Bats=3, Bees=1, Moths=2  (weights ordered Bees, Moths, Bat)
# netchange_pol <- list(
#   list(val = ci_summary(get_net_change(3, w_pol["FDBat_pol"])), lab = "Bats",  col = cols_pol[1], type = "group"),
#   list(type = "spacer", height = 0.5),
#   list(val = ci_summary(get_net_change(1, w_pol["FDBees"])),    lab = "Bees",  col = cols_pol[2], type = "group"),
#   list(type = "spacer", height = 0.5),
#   list(val = ci_summary(get_net_change(2, w_pol["FDMoths"])),   lab = "Moths", col = cols_pol[3], type = "group"),
#   list(type = "spacer", height = 0.5),
#   list(val = ci_summary(get_multi_net_change(c(1,2,3), c(w_pol["FDBees"], w_pol["FDMoths"], w_pol["FDBat_pol"]))),
#        lab = "Multi-group", col = cols_pol[4], type = "multi")
# )

netchange_pol <- list(
  list(val = ci_summary(get_net_change(3)), lab = "Bats",  col = cols_pol[1], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = ci_summary(get_net_change(1)), lab = "Bees",  col = cols_pol[2], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = ci_summary(get_net_change(2)), lab = "Moths", col = cols_pol[3], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = ci_summary(get_multi_net_change(
    c(1,2,3), 
    c(w_pol["FDBees"], w_pol["FDMoths"], w_pol["FDBat_pol"]))),
    lab = "Multi-group", col = cols_pol[4], type = "multi")
)

# Seed dispersal: Bats=4, Birds=5, NF=6
# netchange_sd <- list(
#   list(val = ci_summary(get_net_change(4, w_sd["FDBats"])),  lab = "Bats",              col = cols_sd[1], type = "group"),
#   list(type = "spacer", height = 0.5),
#   list(val = ci_summary(get_net_change(5, w_sd["FDBirds"])), lab = "Birds",             col = cols_sd[2], type = "group"),
#   list(type = "spacer", height = 0.5),
#   list(val = ci_summary(get_net_change(6, w_sd["FDNf"])),    lab = "Non-flying\nmammals", col = cols_sd[3], type = "group"),
#   list(type = "spacer", height = 0.5),
#   list(val = ci_summary(get_multi_net_change(c(4,5,6), c(w_sd["FDBats"], w_sd["FDBirds"], w_sd["FDNf"]))),
#        lab = "Multi-group", col = cols_sd[4], type = "multi")
# )

netchange_sd <- list(
  list(val = ci_summary(get_net_change(4)), lab = "Bats",              col = cols_sd[1], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = ci_summary(get_net_change(5)), lab = "Birds",             col = cols_sd[2], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = ci_summary(get_net_change(6)), lab = "Non-flying\nmammals", col = cols_sd[3], type = "group"),
  list(type = "spacer", height = 0.5),
  list(val = ci_summary(get_multi_net_change(
    c(4,5,6), 
    c(w_sd["FDBats"], w_sd["FDBirds"], w_sd["FDNf"]))),
    lab = "Multi-group", col = cols_sd[4], type = "multi")
)

#### Plotting function (matches draw_rectime_plot style, centred on zero) ####

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

#### X limits (shared across both panels, symmetric around 0) ####

all_vals <- c(
  unlist(lapply(netchange_pol, function(x) if (!is.null(x$val)) x$val else NULL)),
  unlist(lapply(netchange_sd,  function(x) if (!is.null(x$val)) x$val else NULL))
)
x_bound <- max(abs(all_vals), na.rm = TRUE) * 1.1
shared_x_lim <- c(-x_bound, x_bound)

x_bound_pol <- max(abs(unlist(lapply(netchange_pol, function(x) if (!is.null(x$val)) x$val else NULL)))) * 1.1
x_bound_sd <- max(abs(unlist(lapply(netchange_sd, function(x) if (!is.null(x$val)) x$val else NULL)))) * 1.1

x_lim_pol <- c(-x_bound_pol, x_bound_pol)
x_lim_sd <- c(-x_bound_sd, x_bound_sd)

#### Figure ####

svg("output/Figures/Figure_netchange.2.svg", width = 10, height = 4)

layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1))

par(oma = c(1, 1, 1, 1), mgp = c(2, 0.7, 0), tcl = -0.3)

par(mar = c(4, 6, 3, 1))
draw_netchange_plot(netchange_pol, "A) Pollination", x_lim_pol)

par(mar = c(4, 6, 3, 1))
draw_netchange_plot(netchange_sd, "B) Seed dispersal", x_lim_sd)

dev.off()

