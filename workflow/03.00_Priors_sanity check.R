#### Priors sanity check

data <- read.csv("data/processed/model_df.csv")
data$type <- factor(ifelse(1:nrow(data) %in% grep("OG", data$Plot_ID), "old", "rec"),
                    levels = c("old", "rec"))

dataSub <- subset(data, select = c(type, RegTime, ConIndex,  
                                   FDBees, FDMoths,  FDBat_pol, 
                                   FDBats, FDBirds, FDNf
)) 

variable_rec <- colnames(dataSub)[-c(1:3)]
group_index <- 14

# ================================
# PRIOR PREDICTIVE FOR T90
# ================================

set.seed(123)

# ---- mimic your structure ----
n_groups <- length(unique(group_index))
n_var    <- length(group_index)

n_sim <- 4000

# Storage
T90_all <- numeric(n_sim)
lambda_all <- numeric(n_sim)

for(s in 1:n_sim){
  
  # ----- Hyperpriors -----
  
  # Draw one random group
  k <- sample(1:n_groups, 1)
  
  mean_alpha <- rt(1, df = 4) * 0.9 - 2.2
  mean_beta  <- rt(1, df = 4) * 0.35
  
  sigma_alpha <- abs(rt(1, df = 4) * 0.6)
  sigma_beta  <- abs(rt(1, df = 4) * 0.3)
  
  # ----- Draw attribute-level -----
  
  alpha <- rnorm(1, mean_alpha, sigma_alpha)
  beta  <- rnorm(1, mean_beta,  sigma_beta)
  
  # ----- Draw realistic connectivity -----
  # This is crucial: sample from the ACTUAL scaled distribution
  conn <- median(rec_df$connectivity)  
  
  # ----- Rate -----
  lambda <- exp(alpha + beta * conn)
  lambda_all[s] <- lambda
  
  # ----- T90 -----
  T90_all[s] <- log(10) / lambda
}

# ================================
# SUMMARIES
# ================================

quantile(T90_all, c(.01,.05,.1,.25,.5,.75,.9,.95,.99))

# Proportion fast
mean(T90_all < 10)
mean(T90_all < 5)
mean(T90_all < 2)

# Proportion very slow
mean(T90_all > 300)
mean(T90_all > 500)

# Plot
hist(T90_all, breaks=60, xlim=c(0,300),
     main="Prior predictive T90",
     xlab="Years to 90% recovery")

# Also inspect lambda
quantile(lambda_all, c(.01,.05,.5,.95,.99))
                 
plot(lambda_all, T90_all, pch=16, col=rgb(0,0,0,.2),
     xlab="lambda", ylab="T90")
abline(h=c(10,50,200), col="red", lty=2)

# Prior predictive distribution of λ

set.seed(1)

n <- 5000

alpha <- rnorm(n, mean = 0, sd = 2)
beta  <- rnorm(n, mean = 0, sd = 2)
conn  <- rnorm(n, mean = 0, sd = 1)  # your scaled connectivity

lambda <- exp(alpha + beta * conn)

quantile(lambda, c(0.01, 0.1, 0.5, 0.9, 0.99))

# Prior predictive recovery curves

t <- seq(0, 300, by = 1)

theta0  <- 0.2
thetainf <- 1

At <- sapply(lambda[1:200], function(l) {
  theta0 + (thetainf - theta0) * (1 - exp(-l * t))
})

matplot(t, At, type = "l", lty = 1, col = rgb(0,0,0,0.1),
        xlab = "Time", ylab = "Recovery")

# Prior predictive distribution of T90 (most important)
T90 <- rep(NA, 2000)

for(i in 1:2000){
  lam <- exp(rnorm(1, 0, 2) + rnorm(1, 0, 2) * rnorm(1))
  T90[i] <- log(10) / lam
}

quantile(T90, c(0.01, 0.1, 0.5, 0.9, 0.99), na.rm = TRUE)


alpha <- rnorm(n, -2, 1)
beta  <- rnorm(n, 0, 0.3)
conn  <- rnorm(n, 0, 1)
lambda <- exp(alpha + beta * conn)


##############

alpha <- rnorm(5000, -2.2, 0.9)
beta  <- rnorm(5000, 0, 0.35)
conn  <- rnorm(5000, 0, 1)

lambda <- exp(alpha + beta * conn)
T90 <- log(10) / lambda

quantile(T90, c(.05,.1,.25,.5,.75,.9,.95))

