##### Covariance analysis #####

set.seed(1212)

#########################################
#### Pollination and Seed dispersal: ####
#########################################

data <- read.csv("data/processed/model_df.csv")

groups_pol <- c("FDBees", "FDMoths", "FDBat_pol")
groups_disp <- c("FDBats", "FDBirds", "FDNf")
all_groups <- c(groups_pol, groups_disp)
full_groups <- c("FDBees", "FDMoths", "FDBats", "FDBirds")

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

boot_cor <- function(data, indices) {
  d <- data[indices, , drop = FALSE]
  cor(d, use = "pairwise.complete.obs")
}

n <- nrow(net_df_all_z)
R <- 10000
prop <- 0.8

boot_mats_all <- replicate(R, {
  idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
  boot_cor(net_df_all_z, idx)
}, simplify = FALSE)

edge_ci <- function(i, j, mats) {
  vals <- sapply(mats, function(m) m[i, j])
  quantile(vals, c(0.025, 0.975), na.rm = TRUE)
}

p <- ncol(boot_mats_all[[1]])
var_names <- colnames(boot_mats_all[[1]])

edge_ci_arr_all <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_all[i, j, ] <- edge_ci(i, j, boot_mats_all)
    edge_ci_arr_all[j, i, ] <- edge_ci(i, j, boot_mats_all)
  }
}
edge_ci_arr_all

edge_sig_all <- matrix(FALSE, p, p,
                   dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_all)
    edge_sig_all[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_all[j, i] <- edge_sig_all[i, j]
  }
}
edge_sig_all

## Partial correlations controlling for Structure:
# Shared responses to connectivity (and time and other stuff..)

df_resids_str <- data %>% select(all_of(all_groups), "StrIndex") %>%
  mutate(across(all_of(all_groups), ~ log(. + 1)))

get_resids_str <- function(var_name, data) {
  # Usamos na.exclude para que o resíduo mantenha o comprimento original do DF
  form <- as.formula(paste(var_name, "~ StrIndex"))
  mod <- lm(form, data = data, na.action = na.exclude)
  return(resid(mod))
}

net_resids_df_str <- as.data.frame(lapply(all_groups, get_resids_str, data = df_resids_str)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_str) <- all_groups

cor_matrix_res_str <- cor(net_resids_df_str, use = "pairwise.complete.obs")
net_cor_all_ctrl_str <- estimateNetwork(cor_matrix_res_str, 
                                         default = "cor",
                                         labels = colnames(net_resids_df_str))

data_res <- df_resids_str
target_vars <- all_groups
control_vars <- c("StrIndex")
n <- nrow(df_resids_str)

boot_mats_all_ctrl_str <- replicate(R, {
  
  # Sorteia 80% das linhas sem reposição
  idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
  d_sub <- data_res[idx, ]
  
  # Extrair resíduos para cada variável alvo dentro deste sorteio
  resids_list <- lapply(target_vars, function(v) {
    form <- as.formula(paste(v, "~", paste(control_vars, collapse = " + ")))
    mod <- lm(form, data = d_sub, na.action = na.exclude)
    
    # Retornamos o resíduo escalonado (z-score) para manter a estabilidade
    return(as.numeric(scale(resid(mod))))
  })
  
  # Montar matriz de resíduos
  resids_df <- as.data.frame(do.call(cbind, resids_list))
  colnames(resids_df) <- target_vars
  
  # Calcular a correlação dos resíduos
  cor(resids_df, use = "pairwise.complete.obs")
  
}, simplify = FALSE)

p <- ncol(boot_mats_all_ctrl_str[[1]])
var_names <- colnames(boot_mats_all_ctrl_str[[1]])

edge_ci_arr_all_ctrl_str <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_all_ctrl_str[i, j, ] <- edge_ci(i, j, boot_mats_all_ctrl_str)
    edge_ci_arr_all_ctrl_str[j, i, ] <- edge_ci(i, j, boot_mats_all_ctrl_str)
  }
}
edge_ci_arr_all_ctrl_str

edge_sig_all_ctrl_str <- matrix(FALSE, p, p,
                                 dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_all_ctrl_str)
    edge_sig_all_ctrl_str[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_all_ctrl_str[j, i] <- edge_sig_all_ctrl_str[i, j]
  }
}
edge_sig_all_ctrl_str


## Partial correlations controlling for Connectivity:
# Shared responses to structure (and time and other stuff..)

df_resids_con <- data %>% select(all_of(all_groups), "ConIndex") %>%
  mutate(across(all_of(all_groups), ~ log(. + 1)))

get_resids_con <- function(var_name, data) {
  # Usamos na.exclude para que o resíduo mantenha o comprimento original do DF
  form <- as.formula(paste(var_name, "~ ConIndex"))
  mod <- lm(form, data = data, na.action = na.exclude)
  return(resid(mod))
}

net_resids_df_con <- as.data.frame(lapply(all_groups, get_resids_con, data = df_resids_con)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_con) <- all_groups

cor_matrix_res_con <- cor(net_resids_df_con, use = "pairwise.complete.obs")
net_cor_all_ctrl_con <- estimateNetwork(cor_matrix_res_con, 
                                        default = "cor",
                                        labels = colnames(net_resids_df_con))

data_res <- df_resids_con
target_vars <- all_groups
control_vars <- c("ConIndex")
n <- nrow(df_resids_con)

boot_mats_all_ctrl_con <- replicate(R, {
  
  # Sorteia 80% das linhas sem reposição
  idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
  d_sub <- data_res[idx, ]
  
  # Extrair resíduos para cada variável alvo dentro deste sorteio
  resids_list <- lapply(target_vars, function(v) {
    form <- as.formula(paste(v, "~", paste(control_vars, collapse = " + ")))
    mod <- lm(form, data = d_sub, na.action = na.exclude)
    
    # Retornamos o resíduo escalonado (z-score) para manter a estabilidade
    return(as.numeric(scale(resid(mod))))
  })
  
  # Montar matriz de resíduos
  resids_df <- as.data.frame(do.call(cbind, resids_list))
  colnames(resids_df) <- target_vars
  
  # Calcular a correlação dos resíduos
  cor(resids_df, use = "pairwise.complete.obs")
  
}, simplify = FALSE)

p <- ncol(boot_mats_all_ctrl_con[[1]])
var_names <- colnames(boot_mats_all_ctrl_con[[1]])

edge_ci_arr_all_ctrl_con <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_all_ctrl_con[i, j, ] <- edge_ci(i, j, boot_mats_all_ctrl_con)
    edge_ci_arr_all_ctrl_con[j, i, ] <- edge_ci(i, j, boot_mats_all_ctrl_con)
  }
}
edge_ci_arr_all_ctrl_con

edge_sig_all_ctrl_con <- matrix(FALSE, p, p,
                                dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_all_ctrl_con)
    edge_sig_all_ctrl_con[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_all_ctrl_con[j, i] <- edge_sig_all_ctrl_con[i, j]
  }
}
edge_sig_all_ctrl_con


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
net_cor_all_ctrl_both <- estimateNetwork(cor_matrix_res_both, 
                                       default = "cor",
                                       labels = colnames(net_resids_df_both))

data_res <- df_resids_both
target_vars <- all_groups
control_vars <- c("StrIndex", "ConIndex")
n <- nrow(df_resids_both)

boot_mats_all_ctrl_both <- replicate(R, {
  
  # Sorteia 80% das linhas sem reposição
  idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
  d_sub <- data_res[idx, ]
  
  # Extrair resíduos para cada variável alvo dentro deste sorteio
  resids_list <- lapply(target_vars, function(v) {
    form <- as.formula(paste(v, "~", paste(control_vars, collapse = " + ")))
    mod <- lm(form, data = d_sub, na.action = na.exclude)
    
    # Retornamos o resíduo escalonado (z-score) para manter a estabilidade
    return(as.numeric(scale(resid(mod))))
  })
  
  # Montar matriz de resíduos
  resids_df <- as.data.frame(do.call(cbind, resids_list))
  colnames(resids_df) <- target_vars
  
  # Calcular a correlação dos resíduos
  cor(resids_df, use = "pairwise.complete.obs")
  
}, simplify = FALSE)

p <- ncol(boot_mats_all_ctrl_both[[1]])
var_names <- colnames(boot_mats_all_ctrl_both[[1]])

edge_ci_arr_all_ctrl_both <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_all_ctrl_both[i, j, ] <- edge_ci(i, j, boot_mats_all_ctrl_both)
    edge_ci_arr_all_ctrl_both[j, i, ] <- edge_ci(i, j, boot_mats_all_ctrl_both)
  }
}
edge_ci_arr_all_ctrl_both

edge_sig_all_ctrl_both <- matrix(FALSE, p, p,
                       dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_all_ctrl_both)
    edge_sig_all_ctrl_both[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_all_ctrl_both[j, i] <- edge_sig_all_ctrl_both[i, j]
  }
}
edge_sig_all_ctrl_both


###############################
############## Plants #########
###############################

plant_ab <- c("AbSeeds", "AbSdlng")
plant_rich <- c("RichSeeds", "RichSdlng")
plant_shn <- c("ShnSeeds", "ShnSdlng")
plant_fd <- c("FDSeeds", "FDSdlng")

# Abundance

net_df_plants_ab_z <- data %>%
  select(all_of(all_groups), all_of(plant_ab)) %>%
  mutate(across(everything(), ~ scale(log(. + 1))))

cor_matrix_plants_ab <- cor(net_df_plants_ab_z, use = "pairwise.complete.obs", method = "pearson")

net_cor_plants_ab <- estimateNetwork(
  cor_matrix_plants_ab, 
  default = "cor",
  labels = colnames(net_df_plants_ab_z)
)

n <- nrow(net_df_plants_ab_z)
R <- 10000
prop <- 0.8

boot_mats_plants_ab <- replicate(R, {
  idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
  boot_cor(net_df_plants_ab_z, idx)
}, simplify = FALSE)

p <- ncol(boot_mats_plants_ab[[1]])
var_names <- colnames(boot_mats_plants_ab[[1]])

edge_ci_arr_plants_ab <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_plants_ab[i, j, ] <- edge_ci(i, j, boot_mats_plants_ab)
    edge_ci_arr_plants_ab[j, i, ] <- edge_ci(i, j, boot_mats_plants_ab)
  }
}
edge_ci_arr_plants_ab

edge_sig_plants_ab <- matrix(FALSE, p, p,
                       dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_plants_ab)
    edge_sig_plants_ab[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_plants_ab[j, i] <- edge_sig_plants_ab[i, j]
  }
}
edge_sig_plants_ab

# Richness
# Doesnt run. Error:Error in bootnet_correlate(data = data, corMethod = corMethod, corArgs = corArgs,  : 
# Correlation matrix is not positive definite.
# Maybe I can try without the 0s in the future, to see if it runs (0 should be NA)

net_df_plants_rich_z <- data %>%
  select(all_of(full_groups), all_of(plant_rich)) %>%
  mutate(across(everything(), ~ scale(log(. + 1))))

cor_matrix_plants_rich <- cor(net_df_plants_rich_z, use = "pairwise.complete.obs", method = "pearson")

net_cor_plants_rich <- estimateNetwork(
  cor_matrix_plants_rich, 
  default = "cor",
  labels = colnames(net_df_plants_rich_z)
)

n <- nrow(net_df_plants_rich_z)
R <- 10000
prop <- 0.8

boot_mats_plants_rich <- replicate(R, {
  idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
  boot_cor(net_df_plants_rich_z, idx)
}, simplify = FALSE)

p <- ncol(boot_mats_plants_rich[[1]])
var_names <- colnames(boot_mats_plants_rich[[1]])

edge_ci_arr_plants_rich <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_plants_rich[i, j, ] <- edge_ci(i, j, boot_mats_plants_rich)
    edge_ci_arr_plants_rich[j, i, ] <- edge_ci(i, j, boot_mats_plants_rich)
  }
}
edge_ci_arr_plants_rich

edge_sig_plants_rich <- matrix(FALSE, p, p,
                             dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_plants_rich)
    edge_sig_plants_rich[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_plants_rich[j, i] <- edge_sig_plants_rich[i, j]
  }
}
edge_sig_plants_rich

# Shannon
# Doesnt run. Error:Error in bootnet_correlate(data = data, corMethod = corMethod, corArgs = corArgs,  : 
# Correlation matrix is not positive definite.
# Maybe I can try without the 0s in the future, to see if it runs (0 should be NA)

net_df_plants_shn_z <- data %>%
  select(all_of(full_groups), all_of(plant_shn)) %>%
  mutate(across(everything(), ~ scale(log(. + 1))))

cor_matrix_plants_shn <- cor(net_df_plants_shn_z, use = "pairwise.complete.obs", method = "pearson")

# Error in bootnet_correlate(data = data, corMethod = corMethod, corArgs = corArgs,  : 
# Correlation matrix is not positive definite.
# net_cor_plants_shn <- estimateNetwork(
#   cor_matrix_plants_shn, 
#   default = "cor",
#   labels = colnames(net_df_plants_shn_z)
# )
# 
# 
# n <- nrow(net_df_plants_shn_z)
# R <- 10000
# prop <- 0.8
# 
# boot_mats_plants_shn <- replicate(R, {
#   idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
#   boot_cor(net_df_plants_shn_z, idx)
# }, simplify = FALSE)
# 
# p <- ncol(boot_mats_plants_shn[[1]])
# var_names <- colnames(boot_mats_plants_shn[[1]])
# 
# edge_ci_arr_plants_shn <- array(
#   NA,
#   dim = c(p, p, 2),
#   dimnames = list(var_names, var_names, c("low", "high"))
# )
# for (i in 1:(p - 1)) {
#   for (j in (i + 1):p) {
#     edge_ci_arr_plants_shn[i, j, ] <- edge_ci(i, j, boot_mats_plants_shn)
#     edge_ci_arr_plants_shn[j, i, ] <- edge_ci(i, j, boot_mats_plants_shn)
#   }
# }
# edge_ci_arr_plants_shn
# 
# edge_sig_plants_shn <- matrix(FALSE, p, p,
#                              dimnames = list(var_names, var_names))
# for (i in 1:(p - 1)) {
#   for (j in (i + 1):p) {
#     ci <- edge_ci(i, j, boot_mats_plants_shn)
#     edge_sig_plants_shn[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
#     edge_sig_plants_shn[j, i] <- edge_sig_plants_shn[i, j]
#   }
# }
# edge_sig_plants_shn

# Functional Diversity

net_df_plants_fd_z <- data %>%
  select(all_of(full_groups), all_of(plant_fd)) %>%
  mutate(across(everything(), ~ scale(log(. + 1))))

cor_matrix_plants_fd <- cor(net_df_plants_fd_z, use = "pairwise.complete.obs", method = "pearson")

net_cor_plants_fd <- estimateNetwork(
  cor_matrix_plants_fd, 
  default = "cor",
  labels = colnames(net_df_plants_fd_z)
)


n <- nrow(net_df_plants_fd_z)
R <- 10000
prop <- 0.8

boot_mats_plants_fd <- replicate(R, {
  idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
  boot_cor(net_df_plants_fd_z, idx)
}, simplify = FALSE)

p <- ncol(boot_mats_plants_fd[[1]])
var_names <- colnames(boot_mats_plants_fd[[1]])

edge_ci_arr_plants_fd <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_plants_fd[i, j, ] <- edge_ci(i, j, boot_mats_plants_fd)
    edge_ci_arr_plants_fd[j, i, ] <- edge_ci(i, j, boot_mats_plants_fd)
  }
}
edge_ci_arr_plants_fd

edge_sig_plants_fd <- matrix(FALSE, p, p,
                             dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_plants_fd)
    edge_sig_plants_fd[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_plants_fd[j, i] <- edge_sig_plants_fd[i, j]
  }
}
edge_sig_plants_fd


## Partial correlations controlling for Structure:

# Abundance

## Partial correlations controlling for Structure:
# Shared responses to connectivity (and time and other stuff..)

df_resids_plants_ab_str <- data %>% select(all_of(c(all_groups, plant_ab)), "StrIndex") %>%
  mutate(across(all_of(c(all_groups, plant_ab)), ~ log(. + 1)))

net_resids_df_plants_ab_str <- as.data.frame(lapply(c(all_groups, plant_ab), get_resids_str, data = df_resids_plants_ab_str)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_plants_ab_str) <- c(all_groups, plant_ab)

cor_matrix_res_plants_ab_str <- cor(net_resids_df_plants_ab_str, use = "pairwise.complete.obs")
net_cor_plants_ab_ctrl_str <- estimateNetwork(cor_matrix_res_plants_ab_str, 
                                        default = "cor",
                                        labels = colnames(net_resids_df_plants_ab_str))

data_res <- df_resids_plants_ab_str
target_vars <- c(all_groups, plant_ab)
control_vars <- c("StrIndex")
n <- nrow(df_resids_plants_ab_str)

boot_mats_plants_ab_ctrl_str <- replicate(R, {
  
  # Sorteia 80% das linhas sem reposição
  idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
  d_sub <- data_res[idx, ]
  
  # Extrair resíduos para cada variável alvo dentro deste sorteio
  resids_list <- lapply(target_vars, function(v) {
    form <- as.formula(paste(v, "~", paste(control_vars, collapse = " + ")))
    mod <- lm(form, data = d_sub, na.action = na.exclude)
    
    # Retornamos o resíduo escalonado (z-score) para manter a estabilidade
    return(as.numeric(scale(resid(mod))))
  })
  
  # Montar matriz de resíduos
  resids_df <- as.data.frame(do.call(cbind, resids_list))
  colnames(resids_df) <- target_vars
  
  # Calcular a correlação dos resíduos
  cor(resids_df, use = "pairwise.complete.obs")
  
}, simplify = FALSE)

p <- ncol(boot_mats_plants_ab_ctrl_str[[1]])
var_names <- colnames(boot_mats_plants_ab_ctrl_str[[1]])

edge_ci_arr_plants_ab_ctrl_str <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_plants_ab_ctrl_str[i, j, ] <- edge_ci(i, j, boot_mats_plants_ab_ctrl_str)
    edge_ci_arr_plants_ab_ctrl_str[j, i, ] <- edge_ci(i, j, boot_mats_plants_ab_ctrl_str)
  }
}
edge_ci_arr_plants_ab_ctrl_str

edge_sig_plants_ab_ctrl_str <- matrix(FALSE, p, p,
                                dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_plants_ab_ctrl_str)
    edge_sig_plants_ab_ctrl_str[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_plants_ab_ctrl_str[j, i] <- edge_sig_plants_ab_ctrl_str[i, j]
  }
}
edge_sig_plants_ab_ctrl_str

## Partial correlations controlling for Connectivity:
# Shared responses to structure (and time and other stuff..)

df_resids_plants_ab_con <- data %>% select(all_of(c(all_groups, plant_ab)), "ConIndex") %>%
  mutate(across(all_of(c(all_groups, plant_ab)), ~ log(. + 1)))

net_resids_df_plants_ab_con <- as.data.frame(lapply(c(all_groups, plant_ab), get_resids_con, data = df_resids_plants_ab_con)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_plants_ab_con) <- c(all_groups, plant_ab)

cor_matrix_res_plants_ab_con <- cor(net_resids_df_plants_ab_con, use = "pairwise.complete.obs")
net_cor_plants_ab_ctrl_con <- estimateNetwork(cor_matrix_res_plants_ab_con, 
                                              default = "cor",
                                              labels = colnames(net_resids_df_plants_ab_con))

data_res <- df_resids_plants_ab_con
target_vars <- c(all_groups, plant_ab)
control_vars <- c("ConIndex")
n <- nrow(df_resids_plants_ab_con)

boot_mats_plants_ab_ctrl_con <- replicate(R, {
  
  # Sorteia 80% das linhas sem reposição
  idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
  d_sub <- data_res[idx, ]
  
  # Extrair resíduos para cada variável alvo dentro deste sorteio
  resids_list <- lapply(target_vars, function(v) {
    form <- as.formula(paste(v, "~", paste(control_vars, collapse = " + ")))
    mod <- lm(form, data = d_sub, na.action = na.exclude)
    
    # Retornamos o resíduo escalonado (z-score) para manter a estabilidade
    return(as.numeric(scale(resid(mod))))
  })
  
  # Montar matriz de resíduos
  resids_df <- as.data.frame(do.call(cbind, resids_list))
  colnames(resids_df) <- target_vars
  
  # Calcular a correlação dos resíduos
  cor(resids_df, use = "pairwise.complete.obs")
  
}, simplify = FALSE)

p <- ncol(boot_mats_plants_ab_ctrl_con[[1]])
var_names <- colnames(boot_mats_plants_ab_ctrl_con[[1]])

edge_ci_arr_plants_ab_ctrl_con <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_plants_ab_ctrl_con[i, j, ] <- edge_ci(i, j, boot_mats_plants_ab_ctrl_con)
    edge_ci_arr_plants_ab_ctrl_con[j, i, ] <- edge_ci(i, j, boot_mats_plants_ab_ctrl_con)
  }
}
edge_ci_arr_plants_ab_ctrl_con

edge_sig_plants_ab_ctrl_con <- matrix(FALSE, p, p,
                                      dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_plants_ab_ctrl_con)
    edge_sig_plants_ab_ctrl_con[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_plants_ab_ctrl_con[j, i] <- edge_sig_plants_ab_ctrl_con[i, j]
  }
}
edge_sig_plants_ab_ctrl_con

## Partial correlations controlling for Structure and Connectivity:
# Do groups have coordinated recovery once the environment is accounted for?
# shared responses to time and other stuff

df_resids_plants_ab_both <- data %>% select(all_of(c(all_groups, plant_ab)), "StrIndex", "ConIndex") %>%
  mutate(across(all_of(c(all_groups, plant_ab)), ~ log(. + 1)))

net_resids_df_plants_ab_both <- as.data.frame(lapply(c(all_groups, plant_ab), get_resids_both, data = df_resids_plants_ab_both)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_plants_ab_both) <- c(all_groups, plant_ab)

cor_matrix_res_plants_ab_both <- cor(net_resids_df_plants_ab_both, use = "pairwise.complete.obs")
net_cor_plants_ab_ctrl_both <- estimateNetwork(cor_matrix_res_plants_ab_both, 
                                              default = "cor",
                                              labels = colnames(net_resids_df_plants_ab_both))

data_res <- df_resids_plants_ab_both
target_vars <- c(all_groups, plant_ab)
control_vars <- c("StrIndex", "ConIndex")
n <- nrow(df_resids_plants_ab_both)

boot_mats_plants_ab_ctrl_both <- replicate(R, {
  
  # Sorteia 80% das linhas sem reposição
  idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
  d_sub <- data_res[idx, ]
  
  # Extrair resíduos para cada variável alvo dentro deste sorteio
  resids_list <- lapply(target_vars, function(v) {
    form <- as.formula(paste(v, "~", paste(control_vars, collapse = " + ")))
    mod <- lm(form, data = d_sub, na.action = na.exclude)
    
    # Retornamos o resíduo escalonado (z-score) para manter a estabilidade
    return(as.numeric(scale(resid(mod))))
  })
  
  # Montar matriz de resíduos
  resids_df <- as.data.frame(do.call(cbind, resids_list))
  colnames(resids_df) <- target_vars
  
  # Calcular a correlação dos resíduos
  cor(resids_df, use = "pairwise.complete.obs")
  
}, simplify = FALSE)

p <- ncol(boot_mats_plants_ab_ctrl_both[[1]])
var_names <- colnames(boot_mats_plants_ab_ctrl_both[[1]])

edge_ci_arr_plants_ab_ctrl_both <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_plants_ab_ctrl_both[i, j, ] <- edge_ci(i, j, boot_mats_plants_ab_ctrl_both)
    edge_ci_arr_plants_ab_ctrl_both[j, i, ] <- edge_ci(i, j, boot_mats_plants_ab_ctrl_both)
  }
}
edge_ci_arr_plants_ab_ctrl_both

edge_sig_plants_ab_ctrl_both <- matrix(FALSE, p, p,
                                      dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_plants_ab_ctrl_both)
    edge_sig_plants_ab_ctrl_both[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_plants_ab_ctrl_both[j, i] <- edge_sig_plants_ab_ctrl_both[i, j]
  }
}
edge_sig_plants_ab_ctrl_both

# Richness

## Partial correlations controlling for Structure:
# Shared responses to connectivity (and time and other stuff..)

df_resids_plants_rich_str <- data %>% select(all_of(c(full_groups, plant_rich)), "StrIndex") %>%
  mutate(across(all_of(c(full_groups, plant_rich)), ~ log(. + 1)))

net_resids_df_plants_rich_str <- as.data.frame(lapply(c(full_groups, plant_rich), get_resids_str, data = df_resids_plants_rich_str)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_plants_rich_str) <- c(full_groups, plant_rich)

cor_matrix_res_plants_rich_str <- cor(net_resids_df_plants_rich_str, use = "pairwise.complete.obs")
net_cor_plants_rich_ctrl_str <- estimateNetwork(cor_matrix_res_plants_rich_str, 
                                              default = "cor",
                                              labels = colnames(net_resids_df_plants_rich_str))

data_res <- df_resids_plants_rich_str
target_vars <- c(full_groups, plant_rich)
control_vars <- c("StrIndex")
n <- nrow(df_resids_plants_rich_str)

boot_mats_plants_rich_ctrl_str <- replicate(R, {
  
  # Sorteia 80% das linhas sem reposição
  idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
  d_sub <- data_res[idx, ]
  
  # Extrair resíduos para cada variável alvo dentro deste sorteio
  resids_list <- lapply(target_vars, function(v) {
    form <- as.formula(paste(v, "~", paste(control_vars, collapse = " + ")))
    mod <- lm(form, data = d_sub, na.action = na.exclude)
    
    # Retornamos o resíduo escalonado (z-score) para manter a estabilidade
    return(as.numeric(scale(resid(mod))))
  })
  
  # Montar matriz de resíduos
  resids_df <- as.data.frame(do.call(cbind, resids_list))
  colnames(resids_df) <- target_vars
  
  # Calcular a correlação dos resíduos
  cor(resids_df, use = "pairwise.complete.obs")
  
}, simplify = FALSE)

p <- ncol(boot_mats_plants_rich_ctrl_str[[1]])
var_names <- colnames(boot_mats_plants_rich_ctrl_str[[1]])

edge_ci_arr_plants_rich_ctrl_str <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_plants_rich_ctrl_str[i, j, ] <- edge_ci(i, j, boot_mats_plants_rich_ctrl_str)
    edge_ci_arr_plants_rich_ctrl_str[j, i, ] <- edge_ci(i, j, boot_mats_plants_rich_ctrl_str)
  }
}
edge_ci_arr_plants_rich_ctrl_str

edge_sig_plants_rich_ctrl_str <- matrix(FALSE, p, p,
                                      dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_plants_rich_ctrl_str)
    edge_sig_plants_rich_ctrl_str[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_plants_rich_ctrl_str[j, i] <- edge_sig_plants_rich_ctrl_str[i, j]
  }
}
edge_sig_plants_rich_ctrl_str

## Partial correlations controlling for Connectivity:
# Shared responses to structure (and time and other stuff..)
# Doesnt run: Error in bootnet_correlate(data = data, corMethod = corMethod, corArgs = corArgs,  : 
# Correlation matrix is not positive definite.

df_resids_plants_rich_con <- data %>% select(all_of(c(full_groups, plant_rich)), "ConIndex") %>%
  mutate(across(all_of(c(full_groups, plant_rich)), ~ log(. + 1)))

net_resids_df_plants_rich_con <- as.data.frame(lapply(c(full_groups, plant_rich), get_resids_con, data = df_resids_plants_rich_con)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_plants_rich_con) <- c(full_groups, plant_rich)

cor_matrix_res_plants_rich_con <- cor(net_resids_df_plants_rich_con, use = "pairwise.complete.obs")

# Error in bootnet_correlate(data = data, corMethod = corMethod, corArgs = corArgs,  : 
# Correlation matrix is not positive definite.
# net_cor_plants_rich_ctrl_con <- estimateNetwork(cor_matrix_res_plants_rich_con, 
#                                               default = "cor",
#                                               labels = colnames(net_resids_df_plants_rich_con))
# 
# 
# data_res <- df_resids_plants_rich_con
# target_vars <- c(full_groups, plant_rich)
# control_vars <- c("ConIndex")
# n <- nrow(df_resids_plants_rich_con)
# 
# boot_mats_plants_rich_ctrl_con <- replicate(R, {
#   
#   # Sorteia 80% das linhas sem reposição
#   idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
#   d_sub <- data_res[idx, ]
#   
#   # Extrair resíduos para cada variável alvo dentro deste sorteio
#   resids_list <- lapply(target_vars, function(v) {
#     form <- as.formula(paste(v, "~", paste(control_vars, collapse = " + ")))
#     mod <- lm(form, data = d_sub, na.action = na.exclude)
#     
#     # Retornamos o resíduo escalonado (z-score) para manter a estabilidade
#     return(as.numeric(scale(resid(mod))))
#   })
#   
#   # Montar matriz de resíduos
#   resids_df <- as.data.frame(do.call(cbind, resids_list))
#   colnames(resids_df) <- target_vars
#   
#   # Calcular a correlação dos resíduos
#   cor(resids_df, use = "pairwise.complete.obs")
#   
# }, simplify = FALSE)
# 
# p <- ncol(boot_mats_plants_rich_ctrl_con[[1]])
# var_names <- colnames(boot_mats_plants_rich_ctrl_con[[1]])
# 
# edge_ci_arr_plants_rich_ctrl_con <- array(
#   NA,
#   dim = c(p, p, 2),
#   dimnames = list(var_names, var_names, c("low", "high"))
# )
# for (i in 1:(p - 1)) {
#   for (j in (i + 1):p) {
#     edge_ci_arr_plants_rich_ctrl_con[i, j, ] <- edge_ci(i, j, boot_mats_plants_rich_ctrl_con)
#     edge_ci_arr_plants_rich_ctrl_con[j, i, ] <- edge_ci(i, j, boot_mats_plants_rich_ctrl_con)
#   }
# }
# edge_ci_arr_plants_rich_ctrl_con
# 
# edge_sig_plants_rich_ctrl_con <- matrix(FALSE, p, p,
#                                       dimnames = list(var_names, var_names))
# for (i in 1:(p - 1)) {
#   for (j in (i + 1):p) {
#     ci <- edge_ci(i, j, boot_mats_plants_rich_ctrl_con)
#     edge_sig_plants_rich_ctrl_con[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
#     edge_sig_plants_rich_ctrl_con[j, i] <- edge_sig_plants_rich_ctrl_con[i, j]
#   }
# }
# edge_sig_plants_rich_ctrl_con

## Partial correlations controlling for Structure and Connectivity:
# Do groups have coordinated recovery once the environment is accounted for?
# shared responses to time and other stuff

df_resids_plants_rich_both <- data %>% select(all_of(c(full_groups, plant_rich)), "StrIndex", "ConIndex") %>%
  mutate(across(all_of(c(full_groups, plant_rich)), ~ log(. + 1)))

net_resids_df_plants_rich_both <- as.data.frame(lapply(c(full_groups, plant_rich), get_resids_both, data = df_resids_plants_rich_both)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_plants_rich_both) <- c(full_groups, plant_rich)

cor_matrix_res_plants_rich_both <- cor(net_resids_df_plants_rich_both, use = "pairwise.complete.obs")
net_cor_plants_rich_ctrl_both <- estimateNetwork(cor_matrix_res_plants_rich_both, 
                                               default = "cor",
                                               labels = colnames(net_resids_df_plants_rich_both))

data_res <- df_resids_plants_rich_both
target_vars <- c(full_groups, plant_rich)
control_vars <- c("StrIndex", "ConIndex")
n <- nrow(df_resids_plants_rich_both)

boot_mats_plants_rich_ctrl_both <- replicate(R, {
  
  # Sorteia 80% das linhas sem reposição
  idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
  d_sub <- data_res[idx, ]
  
  # Extrair resíduos para cada variável alvo dentro deste sorteio
  resids_list <- lapply(target_vars, function(v) {
    form <- as.formula(paste(v, "~", paste(control_vars, collapse = " + ")))
    mod <- lm(form, data = d_sub, na.action = na.exclude)
    
    # Retornamos o resíduo escalonado (z-score) para manter a estabilidade
    return(as.numeric(scale(resid(mod))))
  })
  
  # Montar matriz de resíduos
  resids_df <- as.data.frame(do.call(cbind, resids_list))
  colnames(resids_df) <- target_vars
  
  # Calcular a correlação dos resíduos
  cor(resids_df, use = "pairwise.complete.obs")
  
}, simplify = FALSE)

p <- ncol(boot_mats_plants_rich_ctrl_both[[1]])
var_names <- colnames(boot_mats_plants_rich_ctrl_both[[1]])

edge_ci_arr_plants_rich_ctrl_both <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_plants_rich_ctrl_both[i, j, ] <- edge_ci(i, j, boot_mats_plants_rich_ctrl_both)
    edge_ci_arr_plants_rich_ctrl_both[j, i, ] <- edge_ci(i, j, boot_mats_plants_rich_ctrl_both)
  }
}
edge_ci_arr_plants_rich_ctrl_both

edge_sig_plants_rich_ctrl_both <- matrix(FALSE, p, p,
                                       dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_plants_rich_ctrl_both)
    edge_sig_plants_rich_ctrl_both[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_plants_rich_ctrl_both[j, i] <- edge_sig_plants_rich_ctrl_both[i, j]
  }
}
edge_sig_plants_rich_ctrl_both

# Shannon

## Partial correlations controlling for Structure:
# Shared responses to connectivity (and time and other stuff..)

df_resids_plants_shn_str <- data %>% select(all_of(c(full_groups, plant_shn)), "StrIndex") %>%
  mutate(across(all_of(c(full_groups, plant_shn)), ~ log(. + 1)))

net_resids_df_plants_shn_str <- as.data.frame(lapply(c(full_groups, plant_shn), get_resids_str, data = df_resids_plants_shn_str)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_plants_shn_str) <- c(full_groups, plant_shn)

cor_matrix_res_plants_shn_str <- cor(net_resids_df_plants_shn_str, use = "pairwise.complete.obs")
net_cor_plants_shn_ctrl_str <- estimateNetwork(cor_matrix_res_plants_shn_str, 
                                                default = "cor",
                                                labels = colnames(net_resids_df_plants_shn_str))

data_res <- df_resids_plants_shn_str
target_vars <- c(full_groups, plant_shn)
control_vars <- c("StrIndex")
n <- nrow(df_resids_plants_shn_str)

boot_mats_plants_shn_ctrl_str <- replicate(R, {
  
  # Sorteia 80% das linhas sem reposição
  idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
  d_sub <- data_res[idx, ]
  
  # Extrair resíduos para cada variável alvo dentro deste sorteio
  resids_list <- lapply(target_vars, function(v) {
    form <- as.formula(paste(v, "~", paste(control_vars, collapse = " + ")))
    mod <- lm(form, data = d_sub, na.action = na.exclude)
    
    # Retornamos o resíduo escalonado (z-score) para manter a estabilidade
    return(as.numeric(scale(resid(mod))))
  })
  
  # Montar matriz de resíduos
  resids_df <- as.data.frame(do.call(cbind, resids_list))
  colnames(resids_df) <- target_vars
  
  # Calcular a correlação dos resíduos
  cor(resids_df, use = "pairwise.complete.obs")
  
}, simplify = FALSE)

p <- ncol(boot_mats_plants_shn_ctrl_str[[1]])
var_names <- colnames(boot_mats_plants_shn_ctrl_str[[1]])

edge_ci_arr_plants_shn_ctrl_str <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_plants_shn_ctrl_str[i, j, ] <- edge_ci(i, j, boot_mats_plants_shn_ctrl_str)
    edge_ci_arr_plants_shn_ctrl_str[j, i, ] <- edge_ci(i, j, boot_mats_plants_shn_ctrl_str)
  }
}
edge_ci_arr_plants_shn_ctrl_str

edge_sig_plants_shn_ctrl_str <- matrix(FALSE, p, p,
                                        dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_plants_shn_ctrl_str)
    edge_sig_plants_shn_ctrl_str[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_plants_shn_ctrl_str[j, i] <- edge_sig_plants_shn_ctrl_str[i, j]
  }
}
edge_sig_plants_shn_ctrl_str

## Partial correlations controlling for Connectivity:
# Shared responses to structure (and time and other stuff..)
# Doesnt run: Error in bootnet_correlate(data = data, corMethod = corMethod, corArgs = corArgs,  : 
# Correlation matrix is not positive definite.

df_resids_plants_shn_con <- data %>% select(all_of(c(full_groups, plant_shn)), "ConIndex") %>%
  mutate(across(all_of(c(full_groups, plant_shn)), ~ log(. + 1)))

net_resids_df_plants_shn_con <- as.data.frame(lapply(c(full_groups, plant_shn), get_resids_con, data = df_resids_plants_shn_con)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_plants_shn_con) <- c(full_groups, plant_shn)

cor_matrix_res_plants_shn_con <- cor(net_resids_df_plants_shn_con, use = "pairwise.complete.obs")

# Error in bootnet_correlate(data = data, corMethod = corMethod, corArgs = corArgs,  : 
# Correlation matrix is not positive definite.
# net_cor_plants_shn_ctrl_con <- estimateNetwork(cor_matrix_res_plants_shn_con, 
#                                                 default = "cor",
#                                                 labels = colnames(net_resids_df_plants_shn_con))
# 
# 
# data_res <- df_resids_plants_shn_con
# target_vars <- c(full_groups, plant_shn)
# control_vars <- c("ConIndex")
# n <- nrow(df_resids_plants_shn_con)
# 
# boot_mats_plants_shn_ctrl_con <- replicate(R, {
#   
#   # Sorteia 80% das linhas sem reposição
#   idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
#   d_sub <- data_res[idx, ]
#   
#   # Extrair resíduos para cada variável alvo dentro deste sorteio
#   resids_list <- lapply(target_vars, function(v) {
#     form <- as.formula(paste(v, "~", paste(control_vars, collapse = " + ")))
#     mod <- lm(form, data = d_sub, na.action = na.exclude)
#     
#     # Retornamos o resíduo escalonado (z-score) para manter a estabilidade
#     return(as.numeric(scale(resid(mod))))
#   })
#   
#   # Montar matriz de resíduos
#   resids_df <- as.data.frame(do.call(cbind, resids_list))
#   colnames(resids_df) <- target_vars
#   
#   # Calcular a correlação dos resíduos
#   cor(resids_df, use = "pairwise.complete.obs")
#   
# }, simplify = FALSE)
# 
# p <- ncol(boot_mats_plants_shn_ctrl_con[[1]])
# var_names <- colnames(boot_mats_plants_shn_ctrl_con[[1]])
# 
# edge_ci_arr_plants_shn_ctrl_con <- array(
#   NA,
#   dim = c(p, p, 2),
#   dimnames = list(var_names, var_names, c("low", "high"))
# )
# for (i in 1:(p - 1)) {
#   for (j in (i + 1):p) {
#     edge_ci_arr_plants_shn_ctrl_con[i, j, ] <- edge_ci(i, j, boot_mats_plants_shn_ctrl_con)
#     edge_ci_arr_plants_shn_ctrl_con[j, i, ] <- edge_ci(i, j, boot_mats_plants_shn_ctrl_con)
#   }
# }
# edge_ci_arr_plants_shn_ctrl_con
# 
# edge_sig_plants_shn_ctrl_con <- matrix(FALSE, p, p,
#                                         dimnames = list(var_names, var_names))
# for (i in 1:(p - 1)) {
#   for (j in (i + 1):p) {
#     ci <- edge_ci(i, j, boot_mats_plants_shn_ctrl_con)
#     edge_sig_plants_shn_ctrl_con[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
#     edge_sig_plants_shn_ctrl_con[j, i] <- edge_sig_plants_shn_ctrl_con[i, j]
#   }
# }
# edge_sig_plants_shn_ctrl_con

## Partial correlations controlling for Structure and Connectivity:
# Do groups have coordinated recovery once the environment is accounted for?
# shared responses to time and other stuff

df_resids_plants_shn_both <- data %>% select(all_of(c(full_groups, plant_shn)), "StrIndex", "ConIndex") %>%
  mutate(across(all_of(c(full_groups, plant_shn)), ~ log(. + 1)))

net_resids_df_plants_shn_both <- as.data.frame(lapply(c(full_groups, plant_shn), get_resids_both, data = df_resids_plants_shn_both)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_plants_shn_both) <- c(full_groups, plant_shn)

cor_matrix_res_plants_shn_both <- cor(net_resids_df_plants_shn_both, use = "pairwise.complete.obs")
net_cor_plants_shn_ctrl_both <- estimateNetwork(cor_matrix_res_plants_shn_both, 
                                                 default = "cor",
                                                 labels = colnames(net_resids_df_plants_shn_both))

data_res <- df_resids_plants_shn_both
target_vars <- c(full_groups, plant_shn)
control_vars <- c("StrIndex", "ConIndex")
n <- nrow(df_resids_plants_shn_both)

boot_mats_plants_shn_ctrl_both <- replicate(R, {
  
  # Sorteia 80% das linhas sem reposição
  idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
  d_sub <- data_res[idx, ]
  
  # Extrair resíduos para cada variável alvo dentro deste sorteio
  resids_list <- lapply(target_vars, function(v) {
    form <- as.formula(paste(v, "~", paste(control_vars, collapse = " + ")))
    mod <- lm(form, data = d_sub, na.action = na.exclude)
    
    # Retornamos o resíduo escalonado (z-score) para manter a estabilidade
    return(as.numeric(scale(resid(mod))))
  })
  
  # Montar matriz de resíduos
  resids_df <- as.data.frame(do.call(cbind, resids_list))
  colnames(resids_df) <- target_vars
  
  # Calcular a correlação dos resíduos
  cor(resids_df, use = "pairwise.complete.obs")
  
}, simplify = FALSE)

p <- ncol(boot_mats_plants_shn_ctrl_both[[1]])
var_names <- colnames(boot_mats_plants_shn_ctrl_both[[1]])

edge_ci_arr_plants_shn_ctrl_both <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_plants_shn_ctrl_both[i, j, ] <- edge_ci(i, j, boot_mats_plants_shn_ctrl_both)
    edge_ci_arr_plants_shn_ctrl_both[j, i, ] <- edge_ci(i, j, boot_mats_plants_shn_ctrl_both)
  }
}
edge_ci_arr_plants_shn_ctrl_both

edge_sig_plants_shn_ctrl_both <- matrix(FALSE, p, p,
                                         dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_plants_shn_ctrl_both)
    edge_sig_plants_shn_ctrl_both[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_plants_shn_ctrl_both[j, i] <- edge_sig_plants_shn_ctrl_both[i, j]
  }
}
edge_sig_plants_shn_ctrl_both

# Functional diversity

## Partial correlations controlling for Structure:
# Shared responses to connectivity (and time and other stuff..)

df_resids_plants_fd_str <- data %>% select(all_of(c(full_groups, plant_fd)), "StrIndex") %>%
  mutate(across(all_of(c(full_groups, plant_fd)), ~ log(. + 1)))

net_resids_df_plants_fd_str <- as.data.frame(lapply(c(full_groups, plant_fd), get_resids_str, data = df_resids_plants_fd_str)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_plants_fd_str) <- c(full_groups, plant_fd)

cor_matrix_res_plants_fd_str <- cor(net_resids_df_plants_fd_str, use = "pairwise.complete.obs")

# Error in bootnet_correlate(data = data, corMethod = corMethod, corArgs = corArgs,  : 
#                              Correlation matrix is not positive definite.
# net_cor_plants_fd_ctrl_str <- estimateNetwork(cor_matrix_res_plants_fd_str, 
#                                                default = "cor",
#                                                labels = colnames(net_resids_df_plants_fd_str))
# 
# data_res <- df_resids_plants_fd_str
# target_vars <- c(full_groups, plant_fd)
# control_vars <- c("StrIndex")
# n <- nrow(df_resids_plants_fd_str)
# 
# boot_mats_plants_fd_ctrl_str <- replicate(R, {
#   
#   # Sorteia 80% das linhas sem reposição
#   idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
#   d_sub <- data_res[idx, ]
#   
#   # Extrair resíduos para cada variável alvo dentro deste sorteio
#   resids_list <- lapply(target_vars, function(v) {
#     form <- as.formula(paste(v, "~", paste(control_vars, collapse = " + ")))
#     mod <- lm(form, data = d_sub, na.action = na.exclude)
#     
#     # Retornamos o resíduo escalonado (z-score) para manter a estabilidade
#     return(as.numeric(scale(resid(mod))))
#   })
#   
#   # Montar matriz de resíduos
#   resids_df <- as.data.frame(do.call(cbind, resids_list))
#   colnames(resids_df) <- target_vars
#   
#   # Calcular a correlação dos resíduos
#   cor(resids_df, use = "pairwise.complete.obs")
#   
# }, simplify = FALSE)
# 
# p <- ncol(boot_mats_plants_fd_ctrl_str[[1]])
# var_names <- colnames(boot_mats_plants_fd_ctrl_str[[1]])
# 
# edge_ci_arr_plants_fd_ctrl_str <- array(
#   NA,
#   dim = c(p, p, 2),
#   dimnames = list(var_names, var_names, c("low", "high"))
# )
# for (i in 1:(p - 1)) {
#   for (j in (i + 1):p) {
#     edge_ci_arr_plants_fd_ctrl_str[i, j, ] <- edge_ci(i, j, boot_mats_plants_fd_ctrl_str)
#     edge_ci_arr_plants_fd_ctrl_str[j, i, ] <- edge_ci(i, j, boot_mats_plants_fd_ctrl_str)
#   }
# }
# edge_ci_arr_plants_fd_ctrl_str
# 
# edge_sig_plants_fd_ctrl_str <- matrix(FALSE, p, p,
#                                        dimnames = list(var_names, var_names))
# for (i in 1:(p - 1)) {
#   for (j in (i + 1):p) {
#     ci <- edge_ci(i, j, boot_mats_plants_fd_ctrl_str)
#     edge_sig_plants_fd_ctrl_str[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
#     edge_sig_plants_fd_ctrl_str[j, i] <- edge_sig_plants_fd_ctrl_str[i, j]
#   }
# }
# edge_sig_plants_fd_ctrl_str

## Partial correlations controlling for Connectivity:
# Shared responses to structure (and time and other stuff..)
# Doesnt run: Error in bootnet_correlate(data = data, corMethod = corMethod, corArgs = corArgs,  : 
# Correlation matrix is not positive definite.

df_resids_plants_fd_con <- data %>% select(all_of(c(full_groups, plant_fd)), "ConIndex") %>%
  mutate(across(all_of(c(full_groups, plant_fd)), ~ log(. + 1)))

net_resids_df_plants_fd_con <- as.data.frame(lapply(c(full_groups, plant_fd), get_resids_con, data = df_resids_plants_fd_con)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_plants_fd_con) <- c(full_groups, plant_fd)

cor_matrix_res_plants_fd_con <- cor(net_resids_df_plants_fd_con, use = "pairwise.complete.obs")

# Error in bootnet_correlate(data = data, corMethod = corMethod, corArgs = corArgs,  : 
#                              Correlation matrix is not positive definite.
# net_cor_plants_fd_ctrl_con <- estimateNetwork(cor_matrix_res_plants_fd_con, 
#                                                default = "cor",
#                                                labels = colnames(net_resids_df_plants_fd_con))
# 
# 
# data_res <- df_resids_plants_fd_con
# target_vars <- c(full_groups, plant_fd)
# control_vars <- c("ConIndex")
# n <- nrow(df_resids_plants_fd_con)
# 
# boot_mats_plants_fd_ctrl_con <- replicate(R, {
#   
#   # Sorteia 80% das linhas sem reposição
#   idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
#   d_sub <- data_res[idx, ]
#   
#   # Extrair resíduos para cada variável alvo dentro deste sorteio
#   resids_list <- lapply(target_vars, function(v) {
#     form <- as.formula(paste(v, "~", paste(control_vars, collapse = " + ")))
#     mod <- lm(form, data = d_sub, na.action = na.exclude)
#     
#     # Retornamos o resíduo escalonado (z-score) para manter a estabilidade
#     return(as.numeric(scale(resid(mod))))
#   })
#   
#   # Montar matriz de resíduos
#   resids_df <- as.data.frame(do.call(cbind, resids_list))
#   colnames(resids_df) <- target_vars
#   
#   # Calcular a correlação dos resíduos
#   cor(resids_df, use = "pairwise.complete.obs")
#   
# }, simplify = FALSE)
# 
# p <- ncol(boot_mats_plants_fd_ctrl_con[[1]])
# var_names <- colnames(boot_mats_plants_fd_ctrl_con[[1]])
# 
# edge_ci_arr_plants_fd_ctrl_con <- array(
#   NA,
#   dim = c(p, p, 2),
#   dimnames = list(var_names, var_names, c("low", "high"))
# )
# for (i in 1:(p - 1)) {
#   for (j in (i + 1):p) {
#     edge_ci_arr_plants_fd_ctrl_con[i, j, ] <- edge_ci(i, j, boot_mats_plants_fd_ctrl_con)
#     edge_ci_arr_plants_fd_ctrl_con[j, i, ] <- edge_ci(i, j, boot_mats_plants_fd_ctrl_con)
#   }
# }
# edge_ci_arr_plants_fd_ctrl_con
# 
# edge_sig_plants_fd_ctrl_con <- matrix(FALSE, p, p,
#                                        dimnames = list(var_names, var_names))
# for (i in 1:(p - 1)) {
#   for (j in (i + 1):p) {
#     ci <- edge_ci(i, j, boot_mats_plants_fd_ctrl_con)
#     edge_sig_plants_fd_ctrl_con[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
#     edge_sig_plants_fd_ctrl_con[j, i] <- edge_sig_plants_fd_ctrl_con[i, j]
#   }
# }
# edge_sig_plants_fd_ctrl_con

## Partial correlations controlling for Structure and Connectivity:
# Do groups have coordinated recovery once the environment is accounted for?
# shared responses to time and other stuff

df_resids_plants_fd_both <- data %>% select(all_of(c(full_groups, plant_fd)), "StrIndex", "ConIndex") %>%
  mutate(across(all_of(c(full_groups, plant_fd)), ~ log(. + 1)))

net_resids_df_plants_fd_both <- as.data.frame(lapply(c(full_groups, plant_fd), get_resids_both, data = df_resids_plants_fd_both)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_plants_fd_both) <- c(full_groups, plant_fd)

cor_matrix_res_plants_fd_both <- cor(net_resids_df_plants_fd_both, use = "pairwise.complete.obs")


# Error in bootnet_correlate(data = data, corMethod = corMethod, corArgs = corArgs,  : 
#                              Correlation matrix is not positive definite.
# net_cor_plants_fd_ctrl_both <- estimateNetwork(cor_matrix_res_plants_fd_both, 
#                                                 default = "cor",
#                                                 labels = colnames(net_resids_df_plants_fd_both))
# 
# data_res <- df_resids_plants_fd_both
# target_vars <- c(full_groups, plant_fd)
# control_vars <- c("StrIndex", "ConIndex")
# n <- nrow(df_resids_plants_fd_both)
# 
# boot_mats_plants_fd_ctrl_both <- replicate(R, {
#   
#   # Sorteia 80% das linhas sem reposição
#   idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
#   d_sub <- data_res[idx, ]
#   
#   # Extrair resíduos para cada variável alvo dentro deste sorteio
#   resids_list <- lapply(target_vars, function(v) {
#     form <- as.formula(paste(v, "~", paste(control_vars, collapse = " + ")))
#     mod <- lm(form, data = d_sub, na.action = na.exclude)
#     
#     # Retornamos o resíduo escalonado (z-score) para manter a estabilidade
#     return(as.numeric(scale(resid(mod))))
#   })
#   
#   # Montar matriz de resíduos
#   resids_df <- as.data.frame(do.call(cbind, resids_list))
#   colnames(resids_df) <- target_vars
#   
#   # Calcular a correlação dos resíduos
#   cor(resids_df, use = "pairwise.complete.obs")
#   
# }, simplify = FALSE)
# 
# p <- ncol(boot_mats_plants_fd_ctrl_both[[1]])
# var_names <- colnames(boot_mats_plants_fd_ctrl_both[[1]])
# 
# edge_ci_arr_plants_fd_ctrl_both <- array(
#   NA,
#   dim = c(p, p, 2),
#   dimnames = list(var_names, var_names, c("low", "high"))
# )
# for (i in 1:(p - 1)) {
#   for (j in (i + 1):p) {
#     edge_ci_arr_plants_fd_ctrl_both[i, j, ] <- edge_ci(i, j, boot_mats_plants_fd_ctrl_both)
#     edge_ci_arr_plants_fd_ctrl_both[j, i, ] <- edge_ci(i, j, boot_mats_plants_fd_ctrl_both)
#   }
# }
# edge_ci_arr_plants_fd_ctrl_both
# 
# edge_sig_plants_fd_ctrl_both <- matrix(FALSE, p, p,
#                                         dimnames = list(var_names, var_names))
# for (i in 1:(p - 1)) {
#   for (j in (i + 1):p) {
#     ci <- edge_ci(i, j, boot_mats_plants_fd_ctrl_both)
#     edge_sig_plants_fd_ctrl_both[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
#     edge_sig_plants_fd_ctrl_both[j, i] <- edge_sig_plants_fd_ctrl_both[i, j]
#   }
# }
# edge_sig_plants_fd_ctrl_both


##### Saving the data:

ci_array_to_df <- function(ci_array) {
  # Extrai os nomes das variáveis
  vars <- dimnames(ci_array)[[1]]
  p <- length(vars)
  
  # Cria todas as combinações de pares (apenas triângulo inferior para evitar repetição)
  indices <- which(lower.tri(matrix(0, p, p)), arr.ind = TRUE)
  
  df <- data.frame(
    Var1  = vars[indices[, 2]],
    Var2  = vars[indices[, 1]],
    Lower = NA,
    Upper = NA
  )
  
  # Preenche com os valores do array
  for(i in 1:nrow(df)) {
    v1 <- df$Var1[i]
    v2 <- df$Var2[i]
    df$Lower[i] <- ci_array[v1, v2, "low"]
    df$Upper[i] <- ci_array[v1, v2, "high"]
  }
  
  return(df)
}

edge_ci_all <- ci_array_to_df(edge_ci_arr_all)
edge_ci_all_ctrl_str <- ci_array_to_df(edge_ci_arr_all_ctrl_str)
edge_ci_all_ctrl_con <- ci_array_to_df(edge_ci_arr_all_ctrl_con)
edge_ci_all_ctrl_both <- ci_array_to_df(edge_ci_arr_all_ctrl_both)

edge_ci_plants_ab <- ci_array_to_df(edge_ci_arr_plants_ab)
edge_ci_plants_ab_ctrl_str <- ci_array_to_df(edge_ci_arr_plants_ab_ctrl_str)
edge_ci_plants_ab_ctrl_con <- ci_array_to_df(edge_ci_arr_plants_ab_ctrl_con)
edge_ci_plants_ab_ctrl_both <- ci_array_to_df(edge_ci_arr_plants_ab_ctrl_both)

edge_ci_plants_rich <- ci_array_to_df(edge_ci_arr_plants_rich)
edge_ci_plants_rich_ctrl_str <- ci_array_to_df(edge_ci_arr_plants_rich_ctrl_str)
# edge_ci_plants_rich_ctrl_con <- ci_array_to_df(edge_ci_arr_plants_rich_ctrl_con)
edge_ci_plants_rich_ctrl_both <- ci_array_to_df(edge_ci_arr_plants_rich_ctrl_both)

# edge_ci_plants_shn <- ci_array_to_df(edge_ci_arr_plants_shn)
edge_ci_plants_shn_ctrl_str <- ci_array_to_df(edge_ci_arr_plants_shn_ctrl_str)
# edge_ci_plants_shn_ctrl_con <- ci_array_to_df(edge_ci_arr_plants_shn_ctrl_con)
edge_ci_plants_shn_ctrl_both <- ci_array_to_df(edge_ci_arr_plants_shn_ctrl_both)

edge_ci_plants_fd <- ci_array_to_df(edge_ci_arr_plants_fd)
# edge_ci_plants_fd_ctrl_str <- ci_array_to_df(edge_ci_arr_plants_fd_ctrl_str)
# edge_ci_plants_fd_ctrl_con <- ci_array_to_df(edge_ci_arr_plants_fd_ctrl_con)
# edge_ci_plants_fd_ctrl_both <- ci_array_to_df(edge_ci_arr_plants_fd_ctrl_both)

results_nets <- list(
  processes = list(
    crude = list(net = net_cor_all, sig = edge_sig_all, ci = edge_ci_all),
    str   = list(net = net_cor_all_ctrl_str, sig = edge_sig_all_ctrl_str, ci = edge_ci_all_ctrl_str),
    con   = list(net = net_cor_all_ctrl_con, sig = edge_sig_all_ctrl_con, ci = edge_ci_all_ctrl_con),
    both  = list(net = net_cor_all_ctrl_both, sig = edge_sig_all_ctrl_both, ci = edge_ci_all_ctrl_both)
  ),
  plants_ab = list(
    crude = list(net = net_cor_plants_ab, sig = edge_sig_plants_ab, ci = edge_ci_plants_ab),
    str   = list(net = net_cor_plants_ab_ctrl_str, sig = edge_sig_plants_ab_ctrl_str, ci = edge_ci_plants_ab_ctrl_str),
    con   = list(net = net_cor_plants_ab_ctrl_con, sig = edge_sig_plants_ab_ctrl_con, ci = edge_ci_plants_ab_ctrl_con),
    both  = list(net = net_cor_plants_ab_ctrl_both, sig = edge_sig_plants_ab_ctrl_both, ci = edge_ci_plants_ab_ctrl_both)
  ),
  plants_rich = list(
    crude = list(net = net_cor_plants_rich, sig = edge_sig_plants_rich, ci = edge_ci_plants_rich),
    str   = list(net = net_cor_plants_rich_ctrl_str, sig = edge_sig_plants_rich_ctrl_str, ci = edge_ci_plants_rich_ctrl_str),
    # con   = list(net = net_cor_plants_rich_ctrl_con, sig = edge_sig_plants_rich_ctrl_con, ci = edge_ci_plants_rich_ctrl_con),
    both  = list(net = net_cor_plants_rich_ctrl_both, sig = edge_sig_plants_rich_ctrl_both, ci = edge_ci_plants_rich_ctrl_both)
  ),
  plants_shn = list(
    # crude = list(net = net_cor_plants_shn, sig = edge_sig_plants_shn, ci = edge_ci_plants_shn),
    str   = list(net = net_cor_plants_shn_ctrl_str, sig = edge_sig_plants_shn_ctrl_str, ci = edge_ci_plants_shn_ctrl_str),
    # con   = list(net = net_cor_plants_shn_ctrl_con, sig = edge_sig_plants_shn_ctrl_con, ci = edge_ci_plants_shn_ctrl_con),
    both  = list(net = net_cor_plants_shn_ctrl_both, sig = edge_sig_plants_shn_ctrl_both, ci = edge_ci_plants_shn_ctrl_both)
  ),
  plants_fd = list(
    crude = list(net = net_cor_plants_fd, sig = edge_sig_plants_fd, ci = edge_ci_plants_fd)#,
    # str   = list(net = net_cor_plants_fd_ctrl_str, sig = edge_sig_plants_fd_ctrl_str, ci = edge_ci_plants_fd_ctrl_str),
    # con   = list(net = net_cor_plants_fd_ctrl_con, sig = edge_sig_plants_fd_ctrl_con, ci = edge_ci_plants_fd_ctrl_con),
    # both  = list(net = net_cor_plants_fd_ctrl_both, sig = edge_sig_plants_fd_ctrl_both, ci = edge_ci_plants_fd_ctrl_both)
  )
)

# saveRDS(results_nets, file = "output/network_results_complete.rds")
