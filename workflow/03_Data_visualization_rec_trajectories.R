
# --- 1. PREPARE RAW DATA ---
# Extract Old Growth data for plotting

var_names <- colnames(dataSub)[-c(1:3)]
raw_og_data <- long %>%
  filter(type == "old") %>%
  mutate(VarName = factor(variable, levels = 1:13, labels = var_names))

# Extract Recovery data
raw_rec_data <- plot_data # From previous steps

# Define the X-axis limit and where to place the OG "Side Panel"
max_time <- max(raw_rec_data$tx, na.rm = TRUE)
gap_width <- max_time * 0.1
og_center <- max_time + gap_width
# Jitter width for OG points
jitter_w  <- max_time * 0.03 

# --- 2. CALCULATE t90 STATISTICS (For the 3rd Panel) ---
# We calculate time to 90% recovery for Low (-1SD), Mean (0), and High (+1SD) Connectivity
conn_levels <- c("Low" = -1, "Mean" = 0, "High" = 1)
t90_list <- list()

# Use a subset of chains for speed
set.seed(123)
sub_chains <- do.call(rbind, samples)[sample(nrow(do.call(rbind, samples)), 500), ]

for(j in 1:9) {
  var_name <- var_names[j]
  
  # Get parameters
  col_t0   <- sub_chains[, grep(paste0("^theta_0\\[", j, "\\]"), colnames(sub_chains))]
  col_tinf <- sub_chains[, grep(paste0("^theta_inf\\[", j, "\\]"), colnames(sub_chains))]
  col_a    <- sub_chains[, grep(paste0("^alpha_con\\[", j, "\\]"), colnames(sub_chains))]
  col_b    <- sub_chains[, grep(paste0("^beta_con\\[", j, "\\]"), colnames(sub_chains))]
  
  for(lvl_name in names(conn_levels)) {
    c_val <- conn_levels[[lvl_name]]
    
    # Calculate Lambda
    lam_vec <- exp(col_a + col_b * c_val)
    
    # Calculate t90: The time when the curve reaches 90% of the gap
    # Equation: 0.9 = 1 - exp(-lambda * t)  ->  0.1 = exp(-lambda * t) 
    # -> ln(0.1) = -lambda * t  ->  t = -ln(0.1) / lambda
    # Note: -ln(0.1) is approx 2.3025
    t90_vec <- -log(0.1) / lam_vec
    
    # Calculate stats
    stats <- quantile(t90_vec, probs = c(0.025, 0.5, 0.975))
    
    t90_list[[paste(j, lvl_name)]] <- data.frame(
      VarName = var_name,
      Connectivity = lvl_name,
      Median = stats[2],
      Lower = stats[1],
      Upper = stats[3]
    )
  }
}
t90_df <- do.call(rbind, t90_list)
t90_df$Connectivity <- factor(t90_df$Connectivity, levels = c("Low", "Mean", "High"))

# --- PLOT A: TRAJECTORIES + OG SIDE PANEL ---

# 1. Prepare Spaghetti & Median Data (Recalculating with the gap logic)
time_seq <- seq(0, max_time, length.out = 100)
# Lista para guardar as iterações (Spaghetti)
full_spag_list <- list()
# For the posterior estimates of OG:
og_est_list <- list()
# Lista para guardar a mediana (Linha Preta)
median_list <- list()

for(j in 1:13) {
  var_name <- var_names[j]
  
  # Extrair parâmetros (usando as sub_chains de 500 amostras que você já tem ou criando novas)
  # Vamos garantir que pegamos 200 para o visual não ficar pesado
  set.seed(123) 
  # Se sub_chains não existir, crie-o novamente a partir de 'samples'
  if(!exists("sub_chains")) sub_chains <- do.call(rbind, samples)[sample(nrow(do.call(rbind, samples)), 200), ]
  
  col_t0   <- sub_chains[, grep(paste0("^theta_0\\[", j, "\\]"), colnames(sub_chains))][1:200]
  col_tinf <- sub_chains[, grep(paste0("^theta_inf\\[", j, "\\]"), colnames(sub_chains))][1:200]
  col_a    <- sub_chains[, grep(paste0("^alpha_con\\[", j, "\\]"), colnames(sub_chains))][1:200]
  col_b    <- sub_chains[, grep(paste0("^beta_con\\[", j, "\\]"), colnames(sub_chains))][1:200]
  
  # Matriz: Linhas = 200 Iterações, Colunas = 100 Pontos de Tempo
  mat_curves <- matrix(NA, nrow=200, ncol=length(time_seq))
  
  for(k in 1:200) {
    # Lambda na conectividade média (0)
    l <- exp(col_a[k] + col_b[k] * 0)
    mat_curves[k, ] <- col_t0[k] + (col_tinf[k] - col_t0[k]) * (1 - exp(-l * time_seq))
  }
  
  # --- 1. Salvar Spaghetti (Formato Longo com coluna 'iter') ---
  iter_df <- as.data.frame(t(mat_curves)) # Transpor para (Tempo x Iter)
  iter_df$tx <- time_seq
  
  # Pivotar para longo
  long_iter <- iter_df %>%
    pivot_longer(cols = -tx, names_to = "iter", values_to = "Y_pred") %>%
    mutate(VarName = var_name)
  
  full_spag_list[[j]] <- long_iter
  
  og_est_list[[j]] <- data.frame(
    VarName = var_names[j],
    tx = og_center,
    Y_med = median(col_tinf),
    Y_lo = quantile(col_tinf, 0.025),
    Y_hi = quantile(col_tinf, 0.975)
  )
  
  # --- 2. Salvar Mediana (Para a linha preta grossa) ---
  med_curve <- apply(mat_curves, 2, median)
  median_list[[j]] <- data.frame(VarName = var_name, tx = time_seq, Y_pred = med_curve)
}


# Consolidar os dataframes
spaghetti_df <- do.call(rbind, full_spag_list) 
median_lines_df <- do.call(rbind, median_list)
og_est_df <- do.call(rbind, og_est_list)

# Garantir que os fatores estão na ordem certa
spaghetti_df$VarName <- factor(spaghetti_df$VarName, levels = var_names)
median_lines_df$VarName <- factor(median_lines_df$VarName, levels = var_names)
og_est_df$VarName <- factor(og_est_df$VarName, levels = var_names)

# 2. CREATE THE MAIN PLOT
p_traj <- ggplot() +
  # --- LEFT SIDE: RECOVERY ---
  
  # Spaghetti (Subset of iterations) - Optional, heavier
  # geom_line(data = spaghetti_df_full, ..., alpha=0.05) + 
  
  # Median Curve
  geom_line(data = median_lines_df, aes(x = tx, y = Y_pred), size = 1) +
  
  # Recovering Data Points (Colored by Connectivity)
  geom_point(data = raw_rec_data, aes(x = tx, y = Y_rec, fill = connectivity), 
             shape = 21, size = 2, alpha = 0.8) +
  
  # --- SEPARATOR LINE ---
  geom_vline(xintercept = max_time + (gap_width/2), linetype = "dotted", color = "grey50") +
  
  # --- RIGHT SIDE: OLD GROWTH ---
  
  # Raw OG Data (Jittered horizontally around og_center)
  geom_jitter(data = raw_og_data, aes(x = og_center, y = value), 
              width = jitter_w, height = 0, 
              shape = 21, fill = "darkgreen", color = "black", alpha = 0.6, size = 2) +
  
  # Model Estimates for OG (PointRange)
  geom_pointrange(data = og_est_df, 
                  aes(x = tx, y = Y_med, ymin = Y_lo, ymax = Y_hi),
                  fill = "gray", shape = 21, size = 0.8, stroke = 1) +
  
  # Labels for the "Fake" X-axis on the right
  annotate("text", x = og_center, y = -Inf, label = "Old\nGrowth", 
           vjust = -0.5, size = 3, fontface = "italic") +
  
  # --- AESTHETICS ---
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, name = "Connectivity") +
  facet_wrap(~VarName, scales = "free_y", ncol = 1, strip.position = "top") + # Stack vertically like base R
  theme_classic() +
  labs(x = "Time (years)", y = "Functional Diversity / Composition", title = "Recovery Trajectories") +
  coord_cartesian(xlim = c(0, og_center + gap_width)) # Clip the view

plot_subset <- function(indices, title_suffix) {
  
  # 1. Identificar variáveis do grupo
  target_vars <- var_names[indices]
  
  # 2. Filtrar os dados para esse grupo
  # AQUI ESTÁ A CORREÇÃO: Filtramos também o spaghetti_df
  sub_spaghetti <- spaghetti_df %>% filter(VarName %in% target_vars)
  
  sub_median    <- median_lines_df %>% filter(VarName %in% target_vars)
  sub_raw_rec   <- raw_rec_data %>% filter(VarName %in% target_vars)
  sub_raw_og    <- raw_og_data %>% filter(VarName %in% target_vars)
  sub_og_est    <- og_est_df %>% filter(VarName %in% target_vars)
  
  # 3. Gerar o Plot
  p <- ggplot() +
    
    # --- A. SPAGHETTI PLOT (AS VÁRIAS LINHAS) ---
    # Alpha bem baixo (0.05 ou 0.1) para criar o efeito de "nuvem"
    geom_line(data = sub_spaghetti, 
              aes(x = tx, y = Y_pred, group = iter), 
              alpha = 0.05, color = "gray50", size = 0.2) +
    
    # --- B. MEDIANA (LINHA GROSSA) ---
    geom_line(data = sub_median, aes(x = tx, y = Y_pred), size = 1.2, color = "black") +
    
    # --- C. PONTOS DE DADOS (RECUPERAÇÃO) ---
    geom_point(data = sub_raw_rec, aes(x = tx, y = Y_rec#, fill = connectivity
                                       ), 
               shape = 21, size = 2, alpha = 0.8) +
    
    # --- D. SEPARADOR E OLD GROWTH ---
    geom_vline(xintercept = max_time + (gap_width/2), linetype = "dotted", color = "grey50") +
    
    geom_jitter(data = sub_raw_og, aes(x = og_center, y = value), 
                width = jitter_w, height = 0, 
                shape = 21, fill = "darkgreen", color = "black", alpha = 0.6, size = 2) +
    
    geom_pointrange(data = sub_og_est, 
                    aes(x = og_center, y = Y_med, ymin = Y_lo, ymax = Y_hi),
                    fill = "gray", shape = 21, size = 0.8, stroke = 1) +
    
    annotate("text", x = og_center, y = -Inf, label = "Old\nGrowth", 
             vjust = -0.5, size = 3, fontface = "italic") +
    
    # --- E. ESTÉTICA ---
    # scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, name = "Connectivity") +
    facet_wrap(~VarName, scales = "free_y", ncol = 1, strip.position = "top") + 
    theme_classic() +
    labs(x = "Time (years)", y = "Functional Diversity") +
    coord_cartesian(xlim = c(0, og_center + gap_width))
  
  return(p)
}

# --- GERANDO OS 3 GRÁFICOS ---

# Grupo 1: Polinizadores (ou o que for 1:3)
plot_group1 <- plot_subset(1:3, "Grupo 1 (Polinizadores)")
print(plot_group1)

# Grupo 2: Dispersores (ou o que for 4:6)
plot_group2 <- plot_subset(4:6, "Grupo 2 (Dispersores)")
print(plot_group2)

# Grupo 4: Estrutura (ou o que for 7:9)
plot_group3 <- plot_subset(7:10, "Grupo 3 (Seedlings)")
print(plot_group3)

# Grupo 4: Estrutura (ou o que for 7:9)
plot_group4 <- plot_subset(11:13, "Grupo 3 (Estrutura)")
print(plot_group4)
