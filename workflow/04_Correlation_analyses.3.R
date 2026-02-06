##### Covariance analysis #####


data <- read.csv("data/processed/model_df.csv")

weights_df <- read.csv("data/processed/weights_df.csv")

w_pol <- c(FDBees=weights_df[1,"prctg"]/sum(weights_df[1:3,"prctg"]), FDMoths=weights_df[2,"prctg"]/sum(weights_df[1:3,"prctg"]), FDBat_pol=weights_df[3,"prctg"]/sum(weights_df[1:3,"prctg"])) 
w_sd <- c(FDBats=weights_df[4,"prctg"]/sum(weights_df[4:6,"prctg"]), FDBirds=weights_df[5,"prctg"]/sum(weights_df[4:6,"prctg"]), FDNf=weights_df[6,"prctg"]/sum(weights_df[4:6,"prctg"]))

data$multi_pol <- apply(data[, names(w_pol)], 1, function(x) {
  sum(x * w_pol, na.rm = TRUE)
})
data$multi_sd <-apply(data[, names(w_sd)], 1, function(x) {
  sum(x * w_sd, na.rm = TRUE)
})

all_groups <- c("multi_pol", "multi_sd")
plant_shn <- c("ShnSeeds", "ShnSdlng")
plant_fd <- c("FDSeeds", "FDSdlng")

# Shannon

net_df_plants_shn_z <- data %>%
  select(all_of(all_groups), all_of(plant_shn)) %>%
  mutate(across(everything(), ~ scale(log(. + 1))))

cor_matrix_plants_shn <- cor(net_df_plants_shn_z, use = "pairwise.complete.obs", method = "pearson")

n <- nrow(net_df_plants_shn_z)
R <- 10000
prop <- 0.8

boot_cor <- function(data, indices) {
  d <- data[indices, , drop = FALSE]
  cor(d, use = "pairwise.complete.obs")
}

boot_mats_plants_shn <- replicate(R, {
  idx <- sample(seq_len(n), size = floor(prop * n), replace = FALSE)
  boot_cor(net_df_plants_shn_z, idx)
}, simplify = FALSE)

p <- ncol(boot_mats_plants_shn[[1]])
var_names <- colnames(boot_mats_plants_shn[[1]])

edge_ci <- function(i, j, mats) {
  vals <- sapply(mats, function(m) m[i, j])
  quantile(vals, c(0.025, 0.975), na.rm = TRUE)
}

edge_ci_arr_plants_shn <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_plants_shn[i, j, ] <- edge_ci(i, j, boot_mats_plants_shn)
    edge_ci_arr_plants_shn[j, i, ] <- edge_ci(i, j, boot_mats_plants_shn)
  }
}
edge_ci_arr_plants_shn

edge_sig_plants_shn <- matrix(FALSE, p, p,
                             dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_plants_shn)
    edge_sig_plants_shn[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_plants_shn[j, i] <- edge_sig_plants_shn[i, j]
  }
}
edge_sig_plants_shn

# Functional Diversity

net_df_plants_fd_z <- data %>%
  select(all_of(all_groups), all_of(plant_fd)) %>%
  mutate(across(everything(), ~ scale(log(. + 1))))

cor_matrix_plants_fd <- cor(net_df_plants_fd_z, use = "pairwise.complete.obs", method = "pearson")


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

# Shannon

## Partial correlations controlling for Structure:
# Shared responses to connectivity (and time and other stuff..)

df_resids_plants_shn_str <- data %>% select(all_of(c(all_groups, plant_shn)), "StrIndex") %>%
  mutate(across(all_of(c(all_groups, plant_shn)), ~ log(. + 1)))

get_resids_str <- function(var_name, data) {
  # Usamos na.exclude para que o resíduo mantenha o comprimento original do DF
  form <- as.formula(paste(var_name, "~ StrIndex"))
  mod <- lm(form, data = data, na.action = na.exclude)
  return(resid(mod))
}

net_resids_df_plants_shn_str <- as.data.frame(lapply(c(all_groups, plant_shn), get_resids_str, data = df_resids_plants_shn_str)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_plants_shn_str) <- c(all_groups, plant_shn)

cor_matrix_res_plants_shn_str <- cor(net_resids_df_plants_shn_str, use = "pairwise.complete.obs")

data_res <- df_resids_plants_shn_str
target_vars <- c(all_groups, plant_shn)
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

df_resids_plants_shn_con <- data %>% select(all_of(c(all_groups, plant_shn)), "ConIndex") %>%
  mutate(across(all_of(c(all_groups, plant_shn)), ~ log(. + 1)))

get_resids_con <- function(var_name, data) {
  # Usamos na.exclude para que o resíduo mantenha o comprimento original do DF
  form <- as.formula(paste(var_name, "~ ConIndex"))
  mod <- lm(form, data = data, na.action = na.exclude)
  return(resid(mod))
}

net_resids_df_plants_shn_con <- as.data.frame(lapply(c(all_groups, plant_shn), get_resids_con, data = df_resids_plants_shn_con)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_plants_shn_con) <- c(all_groups, plant_shn)

cor_matrix_res_plants_shn_con <- cor(net_resids_df_plants_shn_con, use = "pairwise.complete.obs")

data_res <- df_resids_plants_shn_con
target_vars <- c(all_groups, plant_shn)
control_vars <- c("ConIndex")
n <- nrow(df_resids_plants_shn_con)

boot_mats_plants_shn_ctrl_con <- replicate(R, {

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

p <- ncol(boot_mats_plants_shn_ctrl_con[[1]])
var_names <- colnames(boot_mats_plants_shn_ctrl_con[[1]])

edge_ci_arr_plants_shn_ctrl_con <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_plants_shn_ctrl_con[i, j, ] <- edge_ci(i, j, boot_mats_plants_shn_ctrl_con)
    edge_ci_arr_plants_shn_ctrl_con[j, i, ] <- edge_ci(i, j, boot_mats_plants_shn_ctrl_con)
  }
}
edge_ci_arr_plants_shn_ctrl_con

edge_sig_plants_shn_ctrl_con <- matrix(FALSE, p, p,
                                        dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_plants_shn_ctrl_con)
    edge_sig_plants_shn_ctrl_con[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_plants_shn_ctrl_con[j, i] <- edge_sig_plants_shn_ctrl_con[i, j]
  }
}
edge_sig_plants_shn_ctrl_con


# Functional diversity

## Partial correlations controlling for Structure:
# Shared responses to connectivity (and time and other stuff..)

df_resids_plants_fd_str <- data %>% select(all_of(c(all_groups, plant_fd)), "StrIndex") %>%
  mutate(across(all_of(c(all_groups, plant_fd)), ~ log(. + 1)))

net_resids_df_plants_fd_str <- as.data.frame(lapply(c(all_groups, plant_fd), get_resids_str, data = df_resids_plants_fd_str)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_plants_fd_str) <- c(all_groups, plant_fd)

cor_matrix_res_plants_fd_str <- cor(net_resids_df_plants_fd_str, use = "pairwise.complete.obs")

data_res <- df_resids_plants_fd_str
target_vars <- c(all_groups, plant_fd)
control_vars <- c("StrIndex")
n <- nrow(df_resids_plants_fd_str)

boot_mats_plants_fd_ctrl_str <- replicate(R, {

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

p <- ncol(boot_mats_plants_fd_ctrl_str[[1]])
var_names <- colnames(boot_mats_plants_fd_ctrl_str[[1]])

edge_ci_arr_plants_fd_ctrl_str <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_plants_fd_ctrl_str[i, j, ] <- edge_ci(i, j, boot_mats_plants_fd_ctrl_str)
    edge_ci_arr_plants_fd_ctrl_str[j, i, ] <- edge_ci(i, j, boot_mats_plants_fd_ctrl_str)
  }
}
edge_ci_arr_plants_fd_ctrl_str

edge_sig_plants_fd_ctrl_str <- matrix(FALSE, p, p,
                                       dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_plants_fd_ctrl_str)
    edge_sig_plants_fd_ctrl_str[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_plants_fd_ctrl_str[j, i] <- edge_sig_plants_fd_ctrl_str[i, j]
  }
}
edge_sig_plants_fd_ctrl_str

## Partial correlations controlling for Connectivity:
# Shared responses to structure (and time and other stuff..)

df_resids_plants_fd_con <- data %>% select(all_of(c(all_groups, plant_fd)), "ConIndex") %>%
  mutate(across(all_of(c(all_groups, plant_fd)), ~ log(. + 1)))

net_resids_df_plants_fd_con <- as.data.frame(lapply(c(all_groups, plant_fd), get_resids_con, data = df_resids_plants_fd_con)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) # to make sure residuals are in the same scale, so comparable

colnames(net_resids_df_plants_fd_con) <- c(all_groups, plant_fd)

cor_matrix_res_plants_fd_con <- cor(net_resids_df_plants_fd_con, use = "pairwise.complete.obs")

data_res <- df_resids_plants_fd_con
target_vars <- c(all_groups, plant_fd)
control_vars <- c("ConIndex")
n <- nrow(df_resids_plants_fd_con)

boot_mats_plants_fd_ctrl_con <- replicate(R, {

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

p <- ncol(boot_mats_plants_fd_ctrl_con[[1]])
var_names <- colnames(boot_mats_plants_fd_ctrl_con[[1]])

edge_ci_arr_plants_fd_ctrl_con <- array(
  NA,
  dim = c(p, p, 2),
  dimnames = list(var_names, var_names, c("low", "high"))
)
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    edge_ci_arr_plants_fd_ctrl_con[i, j, ] <- edge_ci(i, j, boot_mats_plants_fd_ctrl_con)
    edge_ci_arr_plants_fd_ctrl_con[j, i, ] <- edge_ci(i, j, boot_mats_plants_fd_ctrl_con)
  }
}
edge_ci_arr_plants_fd_ctrl_con

edge_sig_plants_fd_ctrl_con <- matrix(FALSE, p, p,
                                       dimnames = list(var_names, var_names))
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    ci <- edge_ci(i, j, boot_mats_plants_fd_ctrl_con)
    edge_sig_plants_fd_ctrl_con[i, j] <- !(ci[1] <= 0 & ci[2] >= 0)
    edge_sig_plants_fd_ctrl_con[j, i] <- edge_sig_plants_fd_ctrl_con[i, j]
  }
}
edge_sig_plants_fd_ctrl_con


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

edge_ci_plants_shn <- ci_array_to_df(edge_ci_arr_plants_shn)
edge_ci_plants_shn_ctrl_str <- ci_array_to_df(edge_ci_arr_plants_shn_ctrl_str)
edge_ci_plants_shn_ctrl_con <- ci_array_to_df(edge_ci_arr_plants_shn_ctrl_con)

edge_ci_plants_fd <- ci_array_to_df(edge_ci_arr_plants_fd)
edge_ci_plants_fd_ctrl_str <- ci_array_to_df(edge_ci_arr_plants_fd_ctrl_str)
edge_ci_plants_fd_ctrl_con <- ci_array_to_df(edge_ci_arr_plants_fd_ctrl_con)

results_nets <- list(
   plants_shn = list(
    crude = list(net = cor_matrix_plants_shn, sig = edge_sig_plants_shn, ci = edge_ci_plants_shn),
    str   = list(net = cor_matrix_res_plants_shn_str, sig = edge_sig_plants_shn_ctrl_str, ci = edge_ci_plants_shn_ctrl_str),
    con   = list(net = cor_matrix_res_plants_shn_con, sig = edge_sig_plants_shn_ctrl_con, ci = edge_ci_plants_shn_ctrl_con)
  ),
  plants_fd = list(
    crude = list(net = cor_matrix_plants_fd, sig = edge_sig_plants_fd, ci = edge_ci_plants_fd),
    str   = list(net = cor_matrix_res_plants_fd_str, sig = edge_sig_plants_fd_ctrl_str, ci = edge_ci_plants_fd_ctrl_str),
    con   = list(net = cor_matrix_res_plants_fd_con, sig = edge_sig_plants_fd_ctrl_con, ci = edge_ci_plants_fd_ctrl_con)
  )
)

# saveRDS(results_nets, file = "output/network_results_complete.rds")
