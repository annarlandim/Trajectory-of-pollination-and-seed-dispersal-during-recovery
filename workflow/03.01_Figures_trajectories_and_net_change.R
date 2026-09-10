#### Data

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
                                   )) 
dataSub <- dataSub %>%
  mutate(across(
    .cols = -c(1:2),     
    # scaling by dividing the max per column.
    # this is the best for logscale, where negatives are not allowed
    .fns = ~ .x / max(.x, na.rm = TRUE) 
  ))

str(dataSub)

long <- dataSub %>%
  pivot_longer(cols = -c(1:2), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = colnames(dataSub)[-c(1:2)]), variable = as.integer(variable)) %>%
  filter(!is.na(value))

### Figure 2: Trajectories

var_names <- colnames(dataSub)[-c(1,2)]


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

# reordering for plotting

w_pol <- w_pol[c("FDBat_pol", "FDBees", "FDMoths")]

# Attention: Pollinators are ordered as 3, 1, 2
names_pol <- c("Bats", "Bees", "Moths")
names_sd <- c("Bats", "Birds", "Non-flying mammals")

svg("output/Figures/Figure_trajectories.svg", width = 8, height = 5)

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

### Figure 3: Net change

# Posterior of net change

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

# X limits 

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


svg("output/Figures/Figure_netchange.svg", width = 10, height = 4)

layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1))

par(oma = c(1, 1, 1, 1), mgp = c(2, 0.7, 0), tcl = -0.3)

par(mar = c(4, 6, 3, 1))
draw_netchange_plot(netchange_pol, "A) Pollination", x_lim_pol)

par(mar = c(4, 6, 3, 1))
draw_netchange_plot(netchange_sd, "B) Seed dispersal", x_lim_sd)

dev.off()