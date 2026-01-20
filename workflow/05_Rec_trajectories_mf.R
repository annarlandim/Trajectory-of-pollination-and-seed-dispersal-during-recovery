#### Data

qt90 <- readRDS("data/processed/t90.rds")
samples <- readRDS("data/processed/model_posteriors.rds")
mcmc_mat <- as.matrix(samples)

data <- read.csv("data/processed/model_df.csv")
data$type <- factor(ifelse(1:nrow(data) %in% grep("OG", data$Plot_ID), "old", "rec"),
                    levels = c("old", "rec"))

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

library(coda)
# Cores (Mantendo a consistência)
cols_poll <- c("#E41A1C", "#984EA3", "black") # Bees, Moths, BatPol
cols_disp <- c("black", "#377EB8", "#FF7F00") # Bats, Birds, Nf
col_multi <- "#A65628"

# Consolidando apenas Polinização e Dispersão para o Gráfico C
summary_list <- list(
  # Bloco: Polinização
  list(val = qt90$pollination$multi, lab = "Pollination\nMultifunctionality", col = col_multi, type = "multi"),
  list(val = qt90$pollination$per_group[,1], lab = "Bees", col = cols_poll[1], type = "group"),
  list(val = qt90$pollination$per_group[,2], lab = "Moths", col = cols_poll[2], type = "group"),
  list(val = qt90$pollination$per_group[,3], lab = "Bats", col = cols_poll[3], type = "group"),
  
  # Espaço visual
  list(type = "spacer"),
  
  # Bloco: Dispersão
  list(val = qt90$dispersal$multi, lab = "Seed Dispersal\nMultifunctionality", col = col_multi, type = "multi"),
  list(val = qt90$dispersal$per_group[,1], lab = "Bats", col = cols_disp[1], type = "group"),
  list(val = qt90$dispersal$per_group[,2], lab = "Birds", col = cols_disp[2], type = "group"),
  list(val = qt90$dispersal$per_group[,3], lab = "Non-flying mammals", col = cols_disp[3], type = "group")
)

# pdf("recovery_interactions_only.pdf", width = 12, height = 5)

# Layout: A (Trajetória Pol), B (Trajetória Disp), C (Sumário T90)
# Usei uma proporção que dá destaque aos processos individuais e ao sumário
layout_mat <- matrix(c(1, 2, 3, 4, 5), nrow = 1)
layout(layout_mat, widths = c(3, 0.8, 3, 0.8, 4.4))

par(mar = c(4, 4, 3, 0.8), oma = c(1, 1, 1, 1), mgp = c(2, 0.7, 0), tcl = -0.3)

plot_trajectory <- function(var_indices, title, colors, data_subset) {
  plot(NULL, xlim = c(0, 40), ylim = c(0, 1.2), 
       xlab = "Time (years)", ylab = "Functional Diversity", main = title, bty = "l")
  
  for(i in seq_along(var_indices)) {
    idx <- var_indices[i]
    # Subamostra para as linhas de fundo
    sub_lines <- sample(1:nrow(mcmc_mat), 50)
    
    for(s in sub_lines) {
      t0 <- mcmc_mat[s, paste0("theta_0[", idx, "]")]
      tinf <- mcmc_mat[s, paste0("theta_inf[", idx, "]")]
      lamb <- mcmc_mat[s, paste0("lambda[", idx, "]")]
      curve(t0 + (tinf - t0) * (1 - exp(-lamb * x)), 
            add = TRUE, col = adjustcolor(colors[i], alpha.f = 0.1), lwd = 0.5)
    }
    
    # Mediana
    t0_m <- median(mcmc_mat[, paste0("theta_0[", idx, "]")])
    tinf_m <- median(mcmc_mat[, paste0("theta_inf[", idx, "]")])
    lamb_m <- median(mcmc_mat[, paste0("lambda[", idx, "]")])
    curve(t0_m + (tinf_m - t0_m) * (1 - exp(-lamb_m * x)), 
          add = TRUE, col = colors[i], lwd = 3)
    
    # Pontos observados (Recuperação)
    points(data_subset$RegTime[data_subset$variable == idx], 
           data_subset$value[data_subset$variable == idx], 
           pch = 21, bg = adjustcolor(colors[i], alpha.f = 0.6), col = "white", cex = 1.2)
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
plot_trajectory(1:3, "A) Pollination", cols_poll, long[long$type == "rec", ])
par(mar = c(4, 0.8, 3, 1))
plot_old_growth(1:3, cols_poll)

# --- PAINEL B: SEED DISPERSAL ---
par(mar = c(4, 4, 3, 0.8))
plot_trajectory(4:6, "B) Seed Dispersal", cols_disp, long[long$type == "rec", ])
par(mar = c(4, 0.8, 3, 1))
plot_old_growth(4:6, cols_disp)

# --- PAINEL C: RECOVERY SUMMARY (T90) ---
par(mar = c(4, 10, 3, 1))
n_items <- length(summary_list)

# Definir limite do eixo X com base no maior valor de T90 observado nesses grupos
x_max <- max(c(qt90$pollination$per_group, qt90$dispersal$per_group, qt90$pollination$multi, qt90$dispersal$multi), na.rm = TRUE) * 1.1

plot(NULL, xlim = c(0, x_max), ylim = c(0.5, n_items + 0.5), 
     xlab = "Years to 90% recovery", ylab = "", 
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

dev.off()
