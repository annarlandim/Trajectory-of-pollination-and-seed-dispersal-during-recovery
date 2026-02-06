
cor_nets <- readRDS(file = "output/network_results_complete.rds")
str(cor_nets)

# Raw
net_cor_plants_shn <- cor_nets$plants_shn$crude$net
net_cor_plants_fd <- cor_nets$plants_fd$crude$net

edge_sig_plants_shn <- cor_nets$plants_shn$crude$sig
edge_sig_plants_fd <- cor_nets$plants_fd$crude$sig

edge_ci_plants_shn <- cor_nets$plants_shn$crude$ci
edge_ci_plants_fd <- cor_nets$plants_fd$crude$ci

# Controling Forest Structure
net_cor_plants_shn_ctrl_str <- cor_nets$plants_shn$str$net
net_cor_plants_fd_ctrl_str <- cor_nets$plants_fd$str$net

edge_sig_plants_shn_ctrl_str <- cor_nets$plants_shn$str$sig
edge_sig_plants_fd_ctrl_str <- cor_nets$plants_fd$str$sig

edge_ci_plants_shn_ctrl_str <- cor_nets$plants_shn$str$ci
edge_ci_plants_fd_ctrl_str <- cor_nets$plants_fd$str$ci

# Controling Forest Connectivity
net_cor_plants_shn_ctrl_con <- cor_nets$plants_shn$con$net
net_cor_plants_fd_ctrl_con <- cor_nets$plants_fd$con$net

edge_sig_plants_shn_ctrl_con <- cor_nets$plants_shn$con$sig
edge_sig_plants_fd_ctrl_con <- cor_nets$plants_fd$con$sig

edge_ci_plants_shn_ctrl_con <- cor_nets$plants_shn$con$ci
edge_ci_plants_fd_ctrl_con <- cor_nets$plants_fd$con$ci

# Network layouts:

plot_cor_net <- function(cor_matrix, sig_matrix, layout_coords, mode = "full") {
  plot_mat <- cor_matrix
  diag(plot_mat) <- 0 
  
  par(pty = "s")
  
  if (mode == "robust") {
    plot_mat <- plot_mat * sig_matrix 
  }
  
  qgraph(plot_mat, 
         layout = layout_coords,
         vsize = 30,             # Tamanho do nó
         aspect = TRUE,          # Garante que sejam redondos
         labels = FALSE,
         label.cex = 1.2,        # Tamanho do texto DENTRO do nó
         label.scale = FALSE,    # Impede que o texto mude de tamanho entre redes
         rescale = FALSE,         
         posCol = "#018571",      
         negCol = "#a6611a",      
         cut = 0,                 
         mar = c(4, 4, 4, 4),    # Margens internas de cada rede
         theme = "classic")
}
lay_4nodes <- matrix(c(
  -0.8,  1.0,
  0.8,  1.0,
  0.0,  0.4,
  0.0, -0.5
), ncol = 2, byrow = TRUE)

svg(filename = "cor_networks.svg", width = 9, height = 6, pointsize = 12)
par(mfrow = c(2, 3), oma = c(2, 6, 4, 1), mar = c(1, 1, 1, 1))

# --- ROW 1: SHANNON DIVERSITY ---
# Col 1: Raw
plot_cor_net(net_cor_plants_shn, edge_sig_plants_shn, lay_4nodes, mode = "robust")
mtext("Raw", side = 3, line = 1.2, font = 2, cex = 1.1)
mtext("Shannon\nDiversity", side = 2, line = 2.5, font = 2, las = 0, cex = 1)

# Col 2: Structure
plot_cor_net(net_cor_plants_shn_ctrl_str, edge_sig_plants_shn_ctrl_str, lay_4nodes, mode = "robust")
mtext("Structure-controlled", side = 3, line = 1.2, font = 2, cex = 1.1)

# Col 3: Connectivity
plot_cor_net(net_cor_plants_shn_ctrl_con, edge_sig_plants_shn_ctrl_con, lay_4nodes, mode = "robust")
mtext("Connectivity-controlled", side = 3, line = 1.2, font = 2, cex = 1.1)

# --- ROW 2: FUNCTIONAL DIVERSITY ---
# Col 1: Raw
plot_cor_net(net_cor_plants_fd, edge_sig_plants_fd, lay_4nodes, mode = "robust")
mtext("Functional\nDiversity", side = 2, line = 2.5, font = 2, las = 0, cex = 1)

# Col 2: Structure
plot_cor_net(net_cor_plants_fd_ctrl_str, edge_sig_plants_fd_ctrl_str, lay_4nodes, mode = "robust")

# Col 3: Connectivity
plot_cor_net(net_cor_plants_fd_ctrl_con, edge_sig_plants_fd_ctrl_con, lay_4nodes, mode = "robust")

dev.off()
