
cor_nets <- readRDS(file = "output/network_results_complete.rds")
str(cor_nets)

# Raw
net_cor_all <- cor_nets$processes$crude$net
net_cor_plants_ab <- cor_nets$plants_ab$crude$net
net_cor_plants_rich <- cor_nets$plants_rich$crude$net
net_cor_plants_shn <- cor_nets$plants_shn$crude$net
net_cor_plants_fd <- cor_nets$plants_fd$crude$net

edge_sig_all <- cor_nets$processes$crude$sig
edge_sig_plants_ab <- cor_nets$plants_ab$crude$sig
edge_sig_plants_rich <- cor_nets$plants_rich$crude$sig
edge_sig_plants_shn <- cor_nets$plants_shn$crude$sig
edge_sig_plants_fd <- cor_nets$plants_fd$crude$sig

edge_ci_all <- cor_nets$processes$crude$ci
edge_ci_plants_ab <- cor_nets$plants_ab$crude$ci
edge_ci_plants_rich <- cor_nets$plants_rich$crude$ci
edge_ci_plants_shn <- cor_nets$plants_shn$crude$ci
edge_ci_plants_fd <- cor_nets$plants_fd$crude$ci

# Controling Forest Structure
net_cor_all_ctrl_str <- cor_nets$processes$str$net
net_cor_plants_ab_ctrl_str <- cor_nets$plants_ab$str$net
net_cor_plants_rich_ctrl_str <- cor_nets$plants_rich$str$net
net_cor_plants_shn_ctrl_str <- cor_nets$plants_shn$str$net
net_cor_plants_fd_ctrl_str <- cor_nets$plants_fd$str$net

edge_sig_all_ctrl_str <- cor_nets$processes$str$sig
edge_sig_plants_ab_ctrl_str <- cor_nets$plants_ab$str$sig
edge_sig_plants_rich_ctrl_str <- cor_nets$plants_rich$str$sig
edge_sig_plants_shn_ctrl_str <- cor_nets$plants_shn$str$sig
edge_sig_plants_fd_ctrl_str <- cor_nets$plants_fd$str$sig

edge_ci_all_ctrl_str <- cor_nets$processes$str$ci
edge_ci_plants_ab_ctrl_str <- cor_nets$plants_ab$str$ci
edge_ci_plants_rich_ctrl_str <- cor_nets$plants_rich$str$ci
edge_ci_plants_shn_ctrl_str <- cor_nets$plants_shn$str$ci
edge_ci_plants_fd_ctrl_str <- cor_nets$plants_fd$str$ci

# Controling Forest Connectivity
net_cor_all_ctrl_con <- cor_nets$processes$con$net
net_cor_plants_ab_ctrl_con <- cor_nets$plants_ab$con$net
net_cor_plants_rich_ctrl_con <- cor_nets$plants_rich$con$net
net_cor_plants_shn_ctrl_con <- cor_nets$plants_shn$con$net
net_cor_plants_fd_ctrl_con <- cor_nets$plants_fd$con$net

edge_sig_all_ctrl_con <- cor_nets$processes$con$sig
edge_sig_plants_ab_ctrl_con <- cor_nets$plants_ab$con$sig
edge_sig_plants_rich_ctrl_con <- cor_nets$plants_rich$con$sig
edge_sig_plants_shn_ctrl_con <- cor_nets$plants_shn$con$sig
edge_sig_plants_fd_ctrl_con <- cor_nets$plants_fd$con$sig

edge_ci_all_ctrl_con <- cor_nets$processes$con$ci
edge_ci_plants_ab_ctrl_con <- cor_nets$plants_ab$con$ci
edge_ci_plants_rich_ctrl_con <- cor_nets$plants_rich$con$ci
edge_ci_plants_shn_ctrl_con <- cor_nets$plants_shn$con$ci
edge_ci_plants_fd_ctrl_con <- cor_nets$plants_fd$con$ci

# Controling Both
net_cor_all_ctrl_both <- cor_nets$processes$both$net
net_cor_plants_ab_ctrl_both <- cor_nets$plants_ab$both$net
net_cor_plants_rich_ctrl_both <- cor_nets$plants_rich$both$net
net_cor_plants_shn_ctrl_both <- cor_nets$plants_shn$both$net
net_cor_plants_fd_ctrl_both <- cor_nets$plants_fd$both$net

edge_sig_all_ctrl_both <- cor_nets$processes$both$sig
edge_sig_plants_ab_ctrl_both <- cor_nets$plants_ab$both$sig
edge_sig_plants_rich_ctrl_both <- cor_nets$plants_rich$both$sig
edge_sig_plants_shn_ctrl_both <- cor_nets$plants_shn$both$sig
edge_sig_plants_fd_ctrl_both <- cor_nets$plants_fd$both$sig

edge_ci_all_ctrl_both <- cor_nets$processes$both$ci
edge_ci_plants_ab_ctrl_both <- cor_nets$plants_ab$both$ci
edge_ci_plants_rich_ctrl_both <- cor_nets$plants_rich$both$ci
edge_ci_plants_shn_ctrl_both <- cor_nets$plants_shn$both$ci
edge_ci_plants_fd_ctrl_both <- cor_nets$plants_fd$both$ci

# Network layouts:

animal_names <- c(
  "FDBees", "FDMoths", "FDBat_pol",
  "FDBats", "FDBirds", "FDNf"
)

lay_hex_animals <- matrix(c(
  0,  1,    # FDBees
  0.87,  0.5,   # FDMoths
  0.87, -0.5,   # FDBat_pol
  0, -1,    # FDNf
  -0.87, -0.5,   # FDBirds
  -0.87,  0.5    # FDBats
), ncol = 2, byrow = TRUE)

rownames(lay_hex_animals) <- animal_names

plant_names <- c(animal_names, "Seeds", "Seedlings")

lay_hex_plants <- matrix(c(
  0,  1,        # FDBees
  0.87,  0.5,   # FDMoths
  0.87, -0.5,   # FDBat_pol
  0, -1,        # FDNf
  -0.87, -0.5,   # FDBirds
  -0.87,  0.5,   # FDBats
  0,  0.3,      # Seeds (centro cima)
  0, -0.3       # Seedlings (centro baixo)
), ncol = 2, byrow = TRUE)

rownames(lay_hex_plants) <- plant_names

plot_cor_net <- function(cor_matrix, sig_matrix, layout_coords, title = "Network", mode = "full") {
  
  # 1. Prepare the matrix
  plot_mat <- cor_matrix
  diag(plot_mat) <- 0 # Remove self-correlations
  
  # 2. Handle "Robust Only" mode
  if (mode == "robust") {
    # Multiply by logical matrix (TRUE=1, FALSE=0) to keep only significant edges
    plot_mat <- plot_mat * sig_matrix 
  }
  
  # 3. Plotting
  qgraph(plot_mat, 
         layout = layout_coords,
         vsize = 10,              # Corrected from nodesize
         labels = colnames(plot_mat),
         rescale = FALSE,         # Keeps nodes at the exact coordinates provided
         
         # Edge Styling
         posCol = "#2b8cbe",      # Blue for positive
         negCol = "#e41a1c",      # Red for negative
         cut = 0,                 # Starts drawing edges from 0 strength
         
         # Aesthetics
         title = title,
         mar = c(5, 5, 5, 5),
         theme = "classic")
}


# Networks:
par(mfrow=c(1,2))

# Overall networks

## a) Pollination and Seed Dispersal

plot_cor_net(net_cor_all, edge_sig_all, lay_hex_animals, 
             title = "Processes (All)", mode = "full")

plot_cor_net(net_cor_all, edge_sig_all, lay_hex_animals, 
             title = "Processes (Robust)", mode = "robust")

## b) Plants

### Abundance

plot_cor_net(net_cor_plants_ab, edge_sig_plants_ab, lay_hex_plants, 
             title = "Plant Abundance (All)", mode = "full")

plot_cor_net(net_cor_plants_ab, edge_sig_plants_ab, lay_hex_plants, 
             title = "Plant Abundance (Robust", mode = "robust")

### Richness

plot_cor_net(net_cor_plants_rich, edge_sig_plants_rich, lay_hex_plants, 
             title = "Plant Richness (All)", mode = "full")

plot_cor_net(net_cor_plants_rich, edge_sig_plants_rich, lay_hex_plants, 
             title = "Plant Richness (Robust)", mode = "robust")

### Shannon

plot_cor_net(net_cor_plants_shn, edge_sig_plants_shn, lay_hex_plants, 
             title = "Plant Shannon (All)", mode = "full")

plot_cor_net(net_cor_plants_shn, edge_sig_plants_shn, lay_hex_plants, 
             title = "Plant Shannon (Robust)", mode = "robust")

### FD

plot_cor_net(net_cor_plants_fd, edge_sig_plants_fd, lay_hex_plants, 
             title = "Plant FD (All)", mode = "full")

plot_cor_net(net_cor_plants_fd, edge_sig_plants_fd, lay_hex_plants, 
             title = "Plant FD (Robust)", mode = "robust")

# Controlling for structure

## a) Pollination and Seed Dispersal

plot_cor_net(net_cor_all_ctrl_str, edge_sig_all_ctrl_str, lay_hex_animals, 
             title = "Processes (All)", mode = "full")

plot_cor_net(net_cor_all_ctrl_str, edge_sig_all_ctrl_str, lay_hex_animals, 
             title = "Processes (Robust)", mode = "robust")

## b) Plants

### Abundance

plot_cor_net(net_cor_plants_ab_ctrl_str, edge_sig_plants_ab_ctrl_str, lay_hex_plants, 
             title = "Plant Abundance (All)", mode = "full")

plot_cor_net(net_cor_plants_ab_ctrl_str, edge_sig_plants_ab_ctrl_str, lay_hex_plants, 
             title = "Plant Abundance (Robust)", mode = "robust")

### Richness

plot_cor_net(net_cor_plants_rich_ctrl_str, edge_sig_plants_rich_ctrl_str, lay_hex_plants, 
             title = "Plant Richness (All)", mode = "full")

plot_cor_net(net_cor_plants_rich_ctrl_str, edge_sig_plants_rich_ctrl_str, lay_hex_plants, 
             title = "Plant Richness (Robust)", mode = "robust")

### Shannon

plot_cor_net(net_cor_plants_shn_ctrl_str, edge_sig_plants_shn_ctrl_str, lay_hex_plants, 
             title = "Plant Shannon (All)", mode = "full")

plot_cor_net(net_cor_plants_shn_ctrl_str, edge_sig_plants_shn_ctrl_str, lay_hex_plants, 
             title = "Plant Shannon (Robust)", mode = "robust")

### FD

plot_cor_net(net_cor_plants_fd_ctrl_str, edge_sig_plants_fd_ctrl_str, lay_hex_plants, 
             title = "Plant FD (All)", mode = "full")

plot_cor_net(net_cor_plants_fd_ctrl_str, edge_sig_plants_fd_ctrl_str, lay_hex_plants, 
             title = "Plant FD (Robust)", mode = "robust")

# Controlling for connectivity

## a) Pollination and Seed Dispersal

plot_cor_net(net_cor_all_ctrl_con, edge_sig_all_ctrl_con, lay_hex_animals, 
             title = "Processes (All)", mode = "full")

plot_cor_net(net_cor_all_ctrl_con, edge_sig_all_ctrl_con, lay_hex_animals, 
             title = "Processes (Robust)", mode = "robust")

## b) Plants

### Abundance

plot_cor_net(net_cor_plants_ab_ctrl_con, edge_sig_plants_ab_ctrl_con, lay_hex_plants, 
             title = "Plant Abundance (All)", mode = "full")

plot_cor_net(net_cor_plants_ab_ctrl_con, edge_sig_plants_ab_ctrl_con, lay_hex_plants, 
             title = "Plant Abundance (Robust)", mode = "robust")

### Richness

plot_cor_net(net_cor_plants_rich_ctrl_con, edge_sig_plants_rich_ctrl_con, lay_hex_plants, 
             title = "Plant Richness (All)", mode = "full")

plot_cor_net(net_cor_plants_rich_ctrl_con, edge_sig_plants_rich_ctrl_con, lay_hex_plants, 
             title = "Plant Richness (Robust)", mode = "robust")

### Shannon

plot_cor_net(net_cor_plants_shn_ctrl_con, edge_sig_plants_shn_ctrl_con, lay_hex_plants, 
             title = "Plant Shannon (All)", mode = "full")

plot_cor_net(net_cor_plants_shn_ctrl_con, edge_sig_plants_shn_ctrl_con, lay_hex_plants, 
             title = "Plant Shannon (Robust)", mode = "robust")

### FD

plot_cor_net(net_cor_plants_fd_ctrl_con, edge_sig_plants_fd_ctrl_con, lay_hex_plants, 
             title = "Plant FD (All)", mode = "full")

plot_cor_net(net_cor_plants_fd_ctrl_con, edge_sig_plants_fd_ctrl_con, lay_hex_plants, 
             title = "Plant FD (Robust)", mode = "robust")

# Controlling for both

## a) Pollination and Seed Dispersal

plot_cor_net(net_cor_all_ctrl_both, edge_sig_all_ctrl_both, lay_hex_animals, 
             title = "Processes (All)", mode = "full")

plot_cor_net(net_cor_all_ctrl_both, edge_sig_all_ctrl_both, lay_hex_animals, 
             title = "Processes (Robust)", mode = "robust")

## b) Plants

### Abundance

plot_cor_net(net_cor_plants_ab_ctrl_both, edge_sig_plants_ab_ctrl_both, lay_hex_plants, 
             title = "Plant Abundance (All)", mode = "full")

plot_cor_net(net_cor_plants_ab_ctrl_both, edge_sig_plants_ab_ctrl_both, lay_hex_plants, 
             title = "Plant Abundance (Robust)", mode = "robust")

### Richness

plot_cor_net(net_cor_plants_rich_ctrl_both, edge_sig_plants_rich_ctrl_both, lay_hex_plants, 
             title = "Plant Richness (All)", mode = "full")

plot_cor_net(net_cor_plants_rich_ctrl_both, edge_sig_plants_rich_ctrl_both, lay_hex_plants, 
             title = "Plant Richness (Robust)", mode = "robust")

### Shannon

plot_cor_net(net_cor_plants_shn_ctrl_both, edge_sig_plants_shn_ctrl_both, lay_hex_plants, 
             title = "Plant Shannon (All)", mode = "full")

plot_cor_net(net_cor_plants_shn_ctrl_both, edge_sig_plants_shn_ctrl_both, lay_hex_plants, 
             title = "Plant Shannon (Robust)", mode = "robust")

### FD

plot_cor_net(net_cor_plants_fd_ctrl_both, edge_sig_plants_fd_ctrl_both, lay_hex_plants, 
             title = "Plant FD (All)", mode = "full")

plot_cor_net(net_cor_plants_fd_ctrl_both, edge_sig_plants_fd_ctrl_both, lay_hex_plants, 
             title = "Plant FD (Robust)", mode = "robust")







