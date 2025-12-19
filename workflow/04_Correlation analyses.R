##### Covariance analysis #####





# Importatn: in Poorter, 


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
