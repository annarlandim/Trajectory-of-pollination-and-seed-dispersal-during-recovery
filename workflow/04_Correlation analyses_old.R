##### Covariance analysis #####

set.seed(1212)

## Pollination and Seed dispersal:

data <- read.csv("data/processed/model_df.csv")

groups_pol <- c("FDBees", "FDMoths", "FDBat_pol")
groups_disp <- c("FDBats", "FDBirds", "FDNf")
all_groups <- c(groups_pol, groups_disp)

net_df_all <- data[, all_groups]
nrow(net_df_all) # 62

nrow(na.omit(net_df_all)) # 11

net_df_all_z <- net_df_all %>%
  select(all_of(all_groups)) %>%
  mutate(across(everything(), ~ scale(log(. + 1)))) # although no 0s in the data, +1 estabilizes the analysis because log(0>x>1)is really large, creating outliers

# Using pearson correlation, instead of lasso, because allows NAs
# Sparese partial correlation networks (what Poorter uses) requires more data
# "Which groups tend to increase or decrease together during recovery?"
# shared responses to anything

cor_matrix_all <- cor(net_df_all_z, use = "pairwise.complete.obs", method = "pearson")

net_cor_all <- estimateNetwork(
  cor_matrix_all, 
  default = "cor",
  labels = colnames(net_df_all_z)
)
# to test robustness of correlations:
# this alternative is due to having missing data for some groups, different n per group
# the bootnet function does not work, because too many NAs for Nf and Bat_pol
# "2 nodes produced errors; first error: Maximum number of errors in bootstraps reached"
# I tried with only the 4 fuller groups and had the same issue
# Using the package boot:

boot_cor <- function(data, indices) {
  d <- data[indices, ]
  cor(d, use = "pairwise.complete.obs")
}
# bias gives the difference between the mean value of the simulations and the original network
boot_net_all <- boot(
  data = net_df_all_z,
  statistic = boot_cor,
  R = 5000
  )

# Intendfying links:

edge_pairs_net_all <- t(combn(all_groups,2))

get_edge_stats <- function(v1, v2, names_vec, boot_obj) {

  row_idx <- which(names_vec == v1)
  col_idx <- which(names_vec == v2)
  
  # Converter posição da matriz para o índice do vetor 't' (column-major order)
  # Fórmula: (coluna - 1) * total_de_linhas + linha
  idx <- (col_idx - 1) * length(names_vec) + row_idx
  
  # Extrair valores do bootstrap para este índice
  boot_values <- boot_obj$t[, idx]
  
  # Calcular estatísticas
  ci <- quantile(boot_values, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  return(data.frame(
    Var1 = v1,
    Var2 = v2,
    Estimate = ci[2],      # Mediana do bootstrap
    Lower_2.5 = ci[1],     # Limite inferior
    Upper_97.5 = ci[3],    # Limite superior
    Significant = !between(0, ci[1], ci[3])
  ))
}

final_edges_df_net_all <- map2_dfr(edge_pairs_net_all[,1], edge_pairs_net_all[,2], ~get_edge_stats(.x, .y, all_groups, boot_net_all))

final_edges_df_net_all <- final_edges_df_net_all %>% arrange(desc(abs(Estimate)))
rownames(final_edges_df_net_all) <- NULL
final_edges_df_net_all

layout_fixed <- qgraph::qgraph(cor_matrix_all, DoNotPlot = TRUE)$layout
plot(net_cor_all, layout = layout_fixed)

## Partial correlations controlling for Structure and Connectivity:
# Do groups have coordinated recovery once the environment is accounted for?
# shared responses to time and other stuff

df_resids_both <- data %>% select(all_of(all_groups), "StrIndex", "ConIndex") %>%
  mutate(across(all_of(all_groups), ~ log(. + 1)))

get_resids_both <- function(var_name, data) {
  # Usamos na.exclude para que o resíduo mantenha o comprimento original do DF
  form <- as.formula(paste(var_name, "~ ConIndex + StrIndex"))
  mod <- lm(form, data = data, na.action = na.exclude)
  return(resid(mod))
}

net_resids_df_both <- as.data.frame(lapply(all_groups, get_resids_both, data = df_resids_both)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_both) <- all_groups

cor_matrix_res_both <- cor(net_resids_df_both, use = "pairwise.complete.obs")
net_controlled_both <- estimateNetwork(cor_matrix_res_both, 
                                       default = "cor",
                                       labels = colnames(net_resids_df_both))

# Esta função cria a lógica necessária para o bootstrap de resíduos
boot_ctrld_both <- function(target_vars, control_vars) {
  
  function(data, indices) {
    # 1. Reamostrar dados
    d <- data[indices, ]
    
    # 2. Extrair resíduos para cada variável alvo
    resids_list <- lapply(target_vars, function(v) {
      # Constrói a fórmula dinamicamente: Ex: FDBees ~ StrIndex + ConIndex
      formula_str <- paste(v, "~", paste(control_vars, collapse = " + "))
      
      # Tenta rodar o modelo linear
      mod <- try(lm(as.formula(formula_str), data = d, na.action = na.exclude), silent = TRUE)
      
      # Se der erro (ex: muitos NAs no sorteio), retorna NAs
      if(inherits(mod, "try-error")) return(rep(NA, nrow(d)))
      
      # Retorna resíduos escalonados (Z-score)
      return(as.numeric(scale(resid(mod))))
    })
    
    # 3. Montar matriz de resíduos
    resids_mat <- do.call(cbind, resids_list)
    
    # 4. Calcular correlação (retorna vetor para o boot)
    return(cor(resids_mat, use = "pairwise.complete.obs"))
  }
}

boot_net_res_both <- boot(
  data = df_resids_both,
  statistic = boot_ctrld_both(all_groups, c("ConIndex", "StrIndex")),
  R = 5000
)

edge_pairs_net_res_both <- t(combn(all_groups,2))
final_edges_net_res_both <- map2_dfr(edge_pairs_net_res_both[,1], edge_pairs_net_res_both[,2], ~get_edge_stats(.x, .y, all_groups, boot_net_res_both))

final_edges_net_res_both <- final_edges_net_res_both %>% arrange(desc(abs(Estimate)))

plot(net_controlled_both)

## Partial correlations controlling for Structure:
# Shared responses to connectivity (and time and other stuff..)


## Partial correlations controlling for Connectivity:
# Shared responses to structure (and time and other stuff..)

###############################
############## Plants #########
###############################

plant_ab <- c("AbSeeds", "AbSdlng")
plant_rich <- c("RichSeeds", "RichSdlng")
plant_shn <- c("ShnSeeds", "ShnSdlng")
plant_fd <- c("FDSeeds", "FDSdlng")

all_plants <- c(plant_ab, plant_rich, plant_shn, plant_fd)

net_df_plants_ab_z <- data %>%
  select(all_of(all_groups), all_of(plant_ab)) %>%
  mutate(across(everything(), ~ scale(log(. + 1))))

cor_matrix_plants_ab <- cor(net_df_plants_ab_z, use = "pairwise.complete.obs", method = "pearson")

net_cor_plants_ab <- estimateNetwork(
  cor_matrix_plants_ab, 
  default = "cor",
  labels = colnames(net_df_plants_ab_z)
)

boot_net_plants_ab <- boot(
  data = net_df_plants_ab_z,
  statistic = boot_cor,
  R = 5000
)
# Intendfying links:

edge_pairs_net_plants_ab <- t(combn(c(all_groups, plant_ab),2))
final_edges_net_plants_ab <- map2_dfr(edge_pairs_net_plants_ab[,1], edge_pairs_net_plants_ab[,2], ~get_edge_stats(.x, .y, c(all_groups, plant_ab), boot_net_plants_ab))
final_edges_net_plants_ab <- final_edges_net_plants_ab %>% arrange(desc(abs(Estimate)))

plot(net_cor_plants_ab)

## Partial correlations controlling for Structure and Connectivity:

df_resids_plants_ab_both <- data %>% select(all_of(all_groups), all_of(plant_ab), "StrIndex", "ConIndex") %>%
  mutate(across(c(all_of(all_groups), all_of(plant_ab)), ~ log(. + 1)))

net_resids_df_plants_ab_both <- as.data.frame(lapply(c(all_groups, plant_ab), get_resids_both, data = df_resids_plants_ab_both)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_plants_ab_both) <- c(all_groups, plant_ab)

cor_matrix_res_plants_ab_both <- cor(net_resids_df_plants_ab_both, use = "pairwise.complete.obs")
net_controlled_both_plants_ab <- estimateNetwork(cor_matrix_res_plants_ab_both, 
                                       default = "cor",
                                       labels = colnames(net_resids_df_plants_ab_both))

boot_net_res_both_plants_ab <- boot(
  data = df_resids_plants_ab_both,
  statistic = boot_cor,
  R = 5000
)

edge_pairs_net_res_both_plants_ab <- t(combn(c(all_groups, plant_ab),2))
final_edges_net_res_both_plants_ab <- map2_dfr(edge_pairs_net_res_both_plants_ab[,1], edge_pairs_net_res_both_plants_ab[,2], ~get_edge_stats(.x, .y, c(all_groups, plant_ab), boot_net_res_both_plants_ab))

final_edges_net_res_both_plants_ab <- final_edges_net_res_both_plants_ab %>% arrange(desc(abs(Estimate)))

plot(net_controlled_both_plants_ab)

#####################################

# Ideally, now, I should check for multivariate normality of the data. However,
# since I cannot install the package for that, I will move forward:

net_df_all_z <- net_df_all %>%
  mutate(
    logFDBees_z = scale(log(FDBees)),
    logFDMoths_z = scale(log(FDMoths)),
    logFDBat_pol_z = scale(log(FDBat_pol)),
    logFDBats_z = scale(log(FDBats)),
    logFDBirds_z = scale(log(FDBirds)),
    logFDNf_z = scale(log(FDNf))
) %>%
  select(logFDBees_z, logFDMoths_z, logFDBat_pol_z, logFDBats_z, logFDBirds_z, logFDNf_z)

# pearson correlation:
# includes direct and indirect effects. does not distinguish causality from shared effects
# "these variables recover in the same way"
net_cor_pol <- estimateNetwork(
  net_df_all_z,
  default = "cor"
)
plot(net_cor_pol)
boot_net_cor_pol<- bootnet(net_cor_pol)
summary(boot_net_cor_pol)


# partial correlation:
# each link is a conditional association between two groups, i.e., 
# it controls for all other groups in the network.
# "these variables are directly linked"
# poorter controls for 0.5
net_pcor_pol <- estimateNetwork(
  net_df_all_z,
  default = "EBICglasso",
  tuning = 0.5
)

# In case we decide to only select the groups with more data available:

# sum(is.na(data$FDNf))
# sum(is.na(data$FDBat_pol))

net_df_all.2 <- data[, c("FDBees",
                       "FDMoths",
                       "FDBats",
                       "FDBirds")]
net_df_all.2 <- na.omit(net_df_all.2)
nrow(net_df_all.2) # 49

# cor(net_df_all.2)

net_df_all_z.2 <- net_df_all.2 %>%
  mutate(
    logFDBees_z = scale(log(FDBees)),
    logFDMoths_z = scale(log(FDMoths)),
    logFDBats_z = scale(log(FDBats)),
    logFDBirds_z = scale(log(FDBirds))
  ) %>%
  select(logFDBees_z, logFDMoths_z, logFDBats_z, logFDBirds_z)

# pearson correlation
net_cor_pol.2 <- estimateNetwork(
  net_df_all_z.2,
  default = "cor"
)
plot(net_cor_pol.2)

# partial correlation:
# poorter controls for 0.5
net_pcor_pol.2 <- estimateNetwork(
  net_df_all_z.2,
  default = "EBICglasso",
  tuning = 0.5
)
# An empty network was selected to be the best fitting network: so groups are really affected by
# time, instead of affecting each other.

## Pollination

net_df_pol <- data[, c(
  "FDBees",
  "FDMoths",
  "FDBat_pol"
)]
net_df_pol <- na.omit(net_df_pol)
nrow(net_df_pol) # 24

net_df_sd <- data[, c(
  "FDBats",
  "FDBirds",
  "FDNf"
)]
net_df_sd <- na.omit(net_df_sd)
nrow(net_df_sd) # 25

## Per subgroup

post <- as.matrix(do.call(rbind, samples))  
n_var <- constList$n_var
var_names <- colnames(dataSub)[-c(1:3)]  

get_mat <- function(prefix, n_var, post){
  as.matrix(post[, paste0(prefix, "[", 1:n_var, "]")])
}

theta0   <- get_mat("theta_0",   n_var, post)
thetainf <- get_mat("theta_inf", n_var, post)
alpha    <- get_mat("alpha_con", n_var, post)
beta     <- get_mat("beta_con",  n_var, post)

rel_recovery <- function(t, c, theta0, thetainf, alpha, beta){
  lambda <- exp(alpha + beta * c)
  At  <- theta0 + (thetainf - theta0) * (1 - exp(-lambda * t))
  # proximity to OG: 1 when mu = thetainf
  prox <- 1 - abs(At - thetainf) / thetainf
  prox
}

library(huge)
library(bootnet)

### Netowrks 
times <- c(5, 10, 20, 40)
con_med <- quantile(rec_conn, 0.50, na.rm = TRUE)

# full correlation network


cors_ind <- lapply(times, function(t){
  X <- rel_recovery(t, con_med, theta0, thetainf, alpha, beta)
  colnames(X) <- var_names
  estimate(X, use = "pairwise.complete.obs", method = "pearson")
})
names(cors_ind) <- paste0("t", times)

cors_ind$t40

nets_cor <- lapply(times, function(t){
  X <- rel_recovery(t, con_med, theta0, thetainf, alpha, beta)
  colnames(X) <- var_names
  estimateNetwork(X, default = "cor", corMethod = "cor")
})
plot(nets_cor[[2]])
cors_ind$t40

X <- rel_recovery(40, con_med, theta0, thetainf, alpha, beta)
colnames(X) <- var_names
sds <- apply(X, 2, sd)
which(sds == 0 | !is.finite(sds))

# partial correlation networks at 5, 10, 20 e 40 years
# accounts the variation explained by other attributes
# shows independent, causal links between attributes


# Based on Poorter et al. 2021
nets_ind <- lapply(times, function(t){
  X <- rel_recovery(t, con_med, theta0, thetainf, alpha, beta)  # draws x vars
  colnames(X) <- var_names
  
  #Nonparanormal data transformation because EBICglasso assumes multivariate normality
  X_npn <- huge.npn(X, npn.func = "truncation")
  # Partial-correlation network using graphical Lasso
  estimateNetwork(X_npn, default = "EBICglasso", tuning = 0.5)
})

names(nets_ind) <- paste0("t", times)
plot(nets_ind$t40) 

### Per group

groups <- list(
  Pollinators     = c("FDBees","FDMoths","FDBat_pol"),
  SeedDispersers  = c("FDBats","FDBirds","FDNf"),
  ForestStructure = c("AGB","VerticalVH","MaxTH"),
  Seedlings       = c("FDSdlng","AbSdlng","RichSdlng","ShnSdlng")
)

nets_groups <- lapply(times, function(t){
  X <- rel_recovery(t, con_med, theta0, thetainf, alpha, beta)
  colnames(X) <- var_names
  
  # agrega por grupo (média por draw)
  G <- sapply(groups, function(vs) rowMeans(X[, vs, drop = FALSE], na.rm = TRUE))
  
  G_npn <- huge.npn(G, npn.func = "truncation")
  estimateNetwork(G_npn, default = "EBICglasso", tuning = 0.5)
})

names(nets_groups) <- paste0("t", times)
plot(nets_groups$t20)
