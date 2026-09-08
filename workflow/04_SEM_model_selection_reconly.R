set.seed(1212)

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

data_rec <- data %>% filter(Treatment3 != "old-growth forest")

vars_to_scale <- c("ConIndex", "StrIndex", "multi_pol", "multi_sd", "FDSeeds", "FDSdlng")

data_scaled <- data_rec
data_scaled[vars_to_scale] <- lapply(data_rec[vars_to_scale], scale)

complete_rows <- complete.cases(data_scaled[,c("multi_pol", "multi_sd", "ConIndex", "StrIndex", "FDSeeds", "FDSdlng")]) 
complete_data <- data_scaled[complete_rows,c("multi_pol", "multi_sd", "ConIndex", "StrIndex", "FDSeeds", "FDSdlng")]

## Models:

# Full (All paths)
# bf_pol_A   <- bf(multi_pol ~ ConIndex + StrIndex)
# bf_sd_A    <- bf(multi_sd ~ ConIndex + StrIndex)
# bf_seeds_A <- bf(FDSeeds | mi() ~ multi_pol + multi_sd + ConIndex + StrIndex) 
# bf_sdlng_A <- bf(FDSdlng | mi() ~ mi(FDSeeds) + ConIndex + StrIndex)        
# 
# fit_full <- brm(
#   bf_pol_A + bf_sd_A + bf_seeds_A + bf_sdlng_A + set_rescor(FALSE),
#   data = data_scaled, cores = 4, chains = 4, iter = 4000
# )
# 
# summary(fit_full)

# Models:

## Pollination
fit_pol_1 <- brm(
  multi_pol ~ ConIndex,
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000#,
  # file = "output/SEM_models/fit_pol_1_rec"
)
fit_pol_1 <- readRDS("output/SEM_models/fit_pol_1_rec.RDS")

fit_pol_2 <- brm(
  multi_pol ~ StrIndex,
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000#,
 # file = "output/SEM_models/fit_pol_2_rec"
)
fit_pol_2 <- readRDS("output/SEM_models/fit_pol_2_rec.RDS")

fit_pol_3 <- brm(
  multi_pol ~ ConIndex + StrIndex,
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000#,
  # file = "output/SEM_models/fit_pol_3_rec"
)
fit_pol_3 <- readRDS("output/SEM_models/fit_pol_3_rec.RDS")

fit_pol_0 <- brm(
  multi_pol ~ 1,
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000#,
 # file = "output/SEM_models/fit_pol_0_rec"
)
fit_pol_0 <- readRDS("output/SEM_models/fit_pol_0_rec.RDS")

### Compare:

loo_pol_1 <- loo(fit_pol_1, resp = "multipol")
loo_pol_2 <- loo(fit_pol_2, resp = "multipol")
loo_pol_3 <- loo(fit_pol_3, resp = "multipol")
loo_pol_0 <- loo(fit_pol_0, resp = "multipol")

loo_compare(loo_pol_1, loo_pol_2, loo_pol_3, loo_pol_0)

pol_mw <- model_weights(fit_pol_0, fit_pol_1, fit_pol_2, fit_pol_3, resp = "multipol")
pol_mw <- model_weights(fit_pol_0, fit_pol_1, fit_pol_2, fit_pol_3, resp = "multipol", weights = "pseudobma")

posterior_summary(fit_pol_3)

bf_pol <- bf(multi_pol ~ 1)

## Seed dispersal 
# fit_sd_1 <- brm(
#   multi_sd ~ ConIndex,
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sd_1_rec"
# )
fit_sd_1 <- readRDS("output/SEM_models/fit_sd_1_rec.RDS")

# fit_sd_2 <- brm(
#   multi_sd ~ StrIndex,
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sd_2_rec"
# )
fit_sd_2 <- readRDS("output/SEM_models/fit_sd_2_rec.RDS")

# fit_sd_3 <- brm(
#   multi_sd ~ ConIndex + StrIndex,
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sd_3_rec"
# )
fit_sd_3 <- readRDS("output/SEM_models/fit_sd_3_rec.RDS")

# fit_sd_0 <- brm(
#   multi_sd ~ 1,
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sd_0_rec"
# )
fit_sd_0 <- readRDS("output/SEM_models/fit_sd_0_rec.RDS")


### Compare:
loo_sd_1 <- loo(fit_sd_1, resp = "multisd")
loo_sd_2 <- loo(fit_sd_2, resp = "multisd")
loo_sd_3 <- loo(fit_sd_3, resp = "multisd")
loo_sd_0 <- loo(fit_sd_0, resp = "multisd")

loo_compare(loo_sd_1, loo_sd_2, loo_sd_3, loo_sd_0)

bf_sd <- bf(multi_sd ~ StrIndex)

## Seeds

fit_seeds_0 <- brm(
  bf_pol +
  bf_sd +  
  bf(FDSeeds | mi() ~ 1),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_seeds_0_rec"
)
fit_seeds_0 <- readRDS("output/SEM_models/fit_seeds_0_rec.RDS")

fit_seeds_1 <- brm(
  bf_pol +
    bf_sd +  
  bf(FDSeeds | mi() ~ ConIndex),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_seeds_1_rec"
)
fit_seeds_1 <- readRDS("output/SEM_models/fit_seeds_1_rec.RDS")

fit_seeds_2 <- brm(
  bf_pol +
    bf_sd +  
  bf(FDSeeds | mi() ~ StrIndex),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_seeds_2_rec"
)
fit_seeds_2 <- readRDS("output/SEM_models/fit_seeds_2_rec.RDS")

fit_seeds_3 <- brm(
  bf_pol +
    bf_sd +  
  bf(FDSeeds | mi() ~ ConIndex + StrIndex),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_seeds_3_rec"
)
fit_seeds_3 <- readRDS("output/SEM_models/fit_seeds_3_rec.RDS")

fit_seeds_4 <- brm(
  bf_pol +
    bf_sd +  
  bf(FDSeeds | mi() ~ multi_sd),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_seeds_4_rec"
)
fit_seeds_4 <- readRDS("output/SEM_models/fit_seeds_4_rec.RDS")

fit_seeds_5 <- brm(
  bf_pol +
    bf_sd +  
  bf(FDSeeds | mi() ~ multi_pol),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_seeds_5_rec"
)
fit_seeds_5 <- readRDS("output/SEM_models/fit_seeds_5_rec.RDS")

fit_seeds_6 <- brm(
  bf_pol +
    bf_sd +  
  bf(FDSeeds | mi() ~ multi_pol + multi_sd),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_seeds_6_rec"
)
fit_seeds_6 <- readRDS("output/SEM_models/fit_seeds_6_rec.RDS")

fit_seeds_7 <- brm(
  bf_pol +
    bf_sd +  
  bf(FDSeeds | mi() ~ multi_pol + multi_sd + ConIndex),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_seeds_7_rec"
)
fit_seeds_7 <- readRDS("output/SEM_models/fit_seeds_7_rec.RDS")

fit_seeds_8 <- brm(
  bf_pol +
    bf_sd +  
  bf(FDSeeds | mi() ~ multi_pol + multi_sd + StrIndex),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_seeds_8_rec"
)
fit_seeds_8 <- readRDS("output/SEM_models/fit_seeds_8_rec.RDS")

fit_seeds_9 <- brm(
  bf_pol +
    bf_sd +  
  bf(FDSeeds | mi() ~  ConIndex + StrIndex + multi_pol),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_seeds_9_rec"
)
fit_seeds_9 <- readRDS("output/SEM_models/fit_seeds_9_rec.RDS")

fit_seeds_10 <- brm(
  bf_pol +
    bf_sd +  
  bf(FDSeeds | mi() ~ ConIndex + StrIndex + multi_sd),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_seeds_10_rec"
)
fit_seeds_10 <- readRDS("output/SEM_models/fit_seeds_10_rec.RDS")

fit_seeds_11 <- brm(
  bf_pol +
    bf_sd +  
  bf(FDSeeds | mi() ~  ConIndex  + multi_pol),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_seeds_11_rec"
)
fit_seeds_11 <- readRDS("output/SEM_models/fit_seeds_11_rec.RDS")

fit_seeds_12 <- brm(
  bf_pol +
    bf_sd +  
  bf(FDSeeds | mi() ~  ConIndex  + multi_sd),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_seeds_12_rec"
)
fit_seeds_12 <- readRDS("output/SEM_models/fit_seeds_12_rec.RDS")

fit_seeds_13 <- brm(
  bf_pol +
    bf_sd +  
  bf(FDSeeds | mi() ~  StrIndex  + multi_pol),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_seeds_13_rec"
)
fit_seeds_13 <- readRDS("output/SEM_models/fit_seeds_13_rec.RDS")

fit_seeds_14 <- brm(
  bf_pol +
    bf_sd +  
  bf(FDSeeds | mi() ~  StrIndex  + multi_sd),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_seeds_14_rec"
)
fit_seeds_14 <- readRDS("output/SEM_models/fit_seeds_14_rec.RDS")

fit_seeds_15 <- brm(
  bf_pol +
    bf_sd +  
  bf(FDSeeds | mi() ~ ConIndex + StrIndex  + multi_sd + multi_pol),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_seeds_15_rec"
)
fit_seeds_15 <- readRDS("output/SEM_models/fit_seeds_15_rec.RDS")

### Compare:

loo_seeds_0 <- loo(fit_seeds_0, newdata = complete_data, resp = "FDSeeds")
loo_seeds_1 <- loo(fit_seeds_1, newdata = complete_data, resp = "FDSeeds")
loo_seeds_2 <- loo(fit_seeds_2, newdata = complete_data, resp = "FDSeeds")
loo_seeds_3 <- loo(fit_seeds_3, newdata = complete_data, resp = "FDSeeds")
loo_seeds_4 <- loo(fit_seeds_4, newdata = complete_data, resp = "FDSeeds")
loo_seeds_5 <- loo(fit_seeds_5, newdata = complete_data, resp = "FDSeeds")
loo_seeds_6 <- loo(fit_seeds_6, newdata = complete_data, resp = "FDSeeds")
loo_seeds_7 <- loo(fit_seeds_7, newdata = complete_data, resp = "FDSeeds")
loo_seeds_8 <- loo(fit_seeds_8, newdata = complete_data, resp = "FDSeeds")
loo_seeds_9 <- loo(fit_seeds_9, newdata = complete_data, resp = "FDSeeds")
loo_seeds_10 <- loo(fit_seeds_10, newdata = complete_data, resp = "FDSeeds")
loo_seeds_11 <- loo(fit_seeds_11, newdata = complete_data, resp = "FDSeeds")
loo_seeds_12 <- loo(fit_seeds_12, newdata = complete_data, resp = "FDSeeds")
loo_seeds_13 <- loo(fit_seeds_13, newdata = complete_data, resp = "FDSeeds")
loo_seeds_14 <- loo(fit_seeds_14, newdata = complete_data, resp = "FDSeeds")
loo_seeds_15 <- loo(fit_seeds_15, newdata = complete_data, resp = "FDSeeds")

loo_compare(loo_seeds_0, loo_seeds_1, loo_seeds_2, loo_seeds_3, loo_seeds_4, loo_seeds_5, loo_seeds_6, loo_seeds_7, loo_seeds_8, loo_seeds_9, loo_seeds_10, loo_seeds_11, loo_seeds_12, loo_seeds_13, loo_seeds_14, loo_seeds_15)

bf_seeds <-  bf(FDSeeds | mi() ~ multi_sd)

## Seedlings

fit_sdlng_0 <- brm(
  bf_pol +
    bf_sd +  
    bf_seeds +
  bf(FDSdlng | mi() ~ 1),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_sdlng_0_rec"
)

fit_sdlng_1 <- brm(
  bf_pol +
    bf_sd +  
    bf_seeds +
  bf(FDSdlng | mi() ~ ConIndex),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_sdlng_1_rec"
)

fit_sdlng_2 <- brm(
  bf_pol +
    bf_sd +  
    bf_seeds +
  bf(FDSdlng | mi() ~ StrIndex),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_sdlng_2_rec"
)

fit_sdlng_3 <- brm(
  bf_pol +
    bf_sd +  
    bf_seeds +
  bf(FDSdlng | mi() ~ ConIndex + StrIndex),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_sdlng_3_rec"
)

# fit_sdlng_4 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ multi_sd),
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_4_rec"
# )

# fit_sdlng_5 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ multi_pol),
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_5_rec"
# )

# fit_sdlng_6 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ multi_pol + multi_sd),
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_6_rec"
# )

# fit_sdlng_7 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ multi_pol + multi_sd + ConIndex),
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_7_rec"
# )

# fit_sdlng_8 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ multi_pol + multi_sd + StrIndex),
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_8_rec"
# )

# fit_sdlng_9 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~  ConIndex + StrIndex + multi_pol),
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_9_rec"
# )

# fit_sdlng_10 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ ConIndex + StrIndex + multi_sd),
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_10_rec"
# )

# fit_sdlng_11 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~  ConIndex  + multi_pol),
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_11_rec"
# )

# fit_sdlng_12 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~  ConIndex  + multi_sd),
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_12_rec"
# )

# fit_sdlng_13 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~  StrIndex  + multi_pol),
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_13_rec"
# )

# fit_sdlng_14 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~  StrIndex  + multi_sd),
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_14_rec"
# )

# fit_sdlng_15 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ ConIndex + StrIndex  + multi_sd + multi_pol),
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_15_rec"
# )

fit_sdlng_16 <- brm(
  bf_pol +
    bf_sd +  
    bf_seeds +
  bf(FDSdlng | mi() ~ mi(FDSeeds)), # rescor = FALSE by default
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_sdlng_16_rec"
)

fit_sdlng_17 <- brm(
  bf_pol +
    bf_sd +  
    bf_seeds +
  bf (FDSdlng | mi() ~ mi(FDSeeds) + ConIndex), # rescor = FALSE by default
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_sdlng_17_rec"
)

fit_sdlng_18 <- brm(
  bf_pol +
    bf_sd +  
    bf_seeds +
  bf(FDSdlng | mi() ~ mi(FDSeeds) + StrIndex), # rescor = FALSE by default
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_sdlng_18_rec"
)

fit_sdlng_19 <- brm(
  bf_pol +
    bf_sd +  
    bf_seeds +
  bf(FDSdlng | mi() ~ mi(FDSeeds) + StrIndex + ConIndex), # rescor = FALSE by default
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_sdlng_19_rec"
)

# fit_sdlng_20 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ mi(FDSeeds) + multi_pol), # rescor = FALSE by default
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_20_rec"
# )

# fit_sdlng_21 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ mi(FDSeeds) + multi_sd), # rescor = FALSE by default
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_21_rec"
# )

# fit_sdlng_22 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ mi(FDSeeds) + multi_pol + multi_sd), # rescor = FALSE by default
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_22_rec"
# )

# fit_sdlng_23 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ mi(FDSeeds) + multi_pol + ConIndex), # rescor = FALSE by default
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_23_rec"
# )

# fit_sdlng_24 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ mi(FDSeeds) + multi_sd + ConIndex), # rescor = FALSE by default
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_24_rec"
# )

# fit_sdlng_25 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ mi(FDSeeds) + multi_pol + StrIndex), # rescor = FALSE by default
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_25_rec"
# )

# fit_sdlng_26 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ mi(FDSeeds) + multi_sd + StrIndex), # rescor = FALSE by default
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_26_rec"
# )

# fit_sdlng_27 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ mi(FDSeeds) + multi_pol + ConIndex + StrIndex), # rescor = FALSE by default
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_27_rec"
# )

# fit_sdlng_28 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ mi(FDSeeds) + multi_sd + ConIndex + StrIndex), # rescor = FALSE by default
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_28_rec"
# )

# fit_sdlng_29 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ mi(FDSeeds) + multi_pol + multi_sd + ConIndex), # rescor = FALSE by default
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_29_rec"
# )

# fit_sdlng_30 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ mi(FDSeeds) + multi_pol + multi_sd + StrIndex), # rescor = FALSE by default
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_30_rec"
# )

# fit_sdlng_31 <- brm(
#   bf_pol +
#     bf_sd +  
#     bf_seeds +
#   bf(FDSdlng | mi() ~ mi(FDSeeds) + multi_pol + multi_sd + StrIndex + ConIndex), # rescor = FALSE by default
#   data = data_scaled,
#   chains = 4, cores = 4, iter = 4000,
#   file = "output/SEM_models/fit_sdlng_31_rec"
# )

### Compare:
loo_sdlng_0 <- loo(fit_sdlng_0, newdata = complete_data, resp = "FDSdlng")
loo_sdlng_1 <- loo(fit_sdlng_1, newdata = complete_data, resp = "FDSdlng")
loo_sdlng_2 <- loo(fit_sdlng_2, newdata = complete_data, resp = "FDSdlng")
loo_sdlng_3 <- loo(fit_sdlng_3, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_4 <- loo(fit_sdlng_4, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_5 <- loo(fit_sdlng_5, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_6 <- loo(fit_sdlng_6, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_7 <- loo(fit_sdlng_7, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_8 <- loo(fit_sdlng_8, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_9 <- loo(fit_sdlng_9, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_10 <- loo(fit_sdlng_10, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_11 <- loo(fit_sdlng_11, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_12 <- loo(fit_sdlng_12, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_13 <- loo(fit_sdlng_13, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_14 <- loo(fit_sdlng_14, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_15 <- loo(fit_sdlng_15, newdata = complete_data, resp = "FDSdlng")
loo_sdlng_16 <- loo(fit_sdlng_16, newdata = complete_data, resp = "FDSdlng")
loo_sdlng_17 <- loo(fit_sdlng_17, newdata = complete_data, resp = "FDSdlng")
loo_sdlng_18 <- loo(fit_sdlng_18, newdata = complete_data, resp = "FDSdlng")
loo_sdlng_19 <- loo(fit_sdlng_19, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_20 <- loo(fit_sdlng_20, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_21 <- loo(fit_sdlng_21, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_22 <- loo(fit_sdlng_22, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_23 <- loo(fit_sdlng_23, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_24 <- loo(fit_sdlng_24, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_25 <- loo(fit_sdlng_25, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_26 <- loo(fit_sdlng_26, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_27 <- loo(fit_sdlng_27, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_28 <- loo(fit_sdlng_28, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_29 <- loo(fit_sdlng_29, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_30 <- loo(fit_sdlng_30, newdata = complete_data, resp = "FDSdlng")
# loo_sdlng_31 <- loo(fit_sdlng_31, newdata = complete_data, resp = "FDSdlng")

loo_compare(loo_sdlng_0, loo_sdlng_1, loo_sdlng_2, loo_sdlng_3, 
            # loo_sdlng_4, loo_sdlng_5, loo_sdlng_6, loo_sdlng_7, 
            # loo_sdlng_8, loo_sdlng_9, loo_sdlng_10, loo_sdlng_11, loo_sdlng_12, loo_sdlng_13, loo_sdlng_14, loo_sdlng_15,
            loo_sdlng_16, loo_sdlng_17, loo_sdlng_18, loo_sdlng_19
            # loo_sdlng_20, 
            # loo_sdlng_21, loo_sdlng_22, loo_sdlng_23, loo_sdlng_24, loo_sdlng_25, 
            # loo_sdlng_26, 
            # loo_sdlng_27, 
            # loo_sdlng_28, loo_sdlng_29, loo_sdlng_30, loo_sdlng_31
            )

# Final model:

fit_final <- brm(
  bf(multi_pol ~ StrIndex) +
    bf(multi_sd ~ StrIndex) +  
    bf(FDSeeds | mi() ~ multi_sd) +
    bf(FDSdlng | mi() ~ 1),
  data = data_scaled,
  chains = 4, cores = 4, iter = 4000,
  file = "output/SEM_models/fit_final_rec"
)

summary(fit_final)




# Mediated: Str and con only via pol and sd (B)
# bf_pol_B   <- bf(multi_pol ~ ConIndex + StrIndex)
# bf_sd_B    <- bf(multi_sd ~ ConIndex + StrIndex)
# bf_seeds_B <- bf(FDSeeds | mi() ~ multi_pol + multi_sd) # Sem Con e Str diretos
# bf_sdlng_B <- bf(FDSdlng | mi() ~ mi(FDSeeds))         # Só sementes importam para plântulas
# 
# fit_mediated <- brm(
#   bf_pol_B + bf_sd_B + bf_seeds_B + bf_sdlng_B + set_rescor(FALSE),
#   data = data_scaled, cores = 4, chains = 4, iter = 4000,
#   file = "output/SEM_models/fit_mediated"
# )

fit_mediated <- readRDS("output/SEM_models/fit_mediated.RDS")

summary(fit_mediated)

# Ab+Bio: Str and biotic affect seedlings (C)
# bf_pol_C   <- bf(multi_pol ~ ConIndex + StrIndex)
# bf_sd_C    <- bf(multi_sd ~ ConIndex + StrIndex)
# bf_seeds_C <- bf(FDSeeds | mi() ~ multi_pol + multi_sd) 
# bf_sdlng_C <- bf(FDSdlng | mi() ~ mi(FDSeeds) + StrIndex)         
# 
# fit_abbio <- brm(
#   bf_pol_C + bf_sd_C + bf_seeds_C + bf_sdlng_C + set_rescor(FALSE),
#   data = data_scaled, cores = 4, chains = 4, iter = 4000,
#   file = "output/SEM_models/fit_abbio"
# )

fit_abbio <- readRDS("output/SEM_models/fit_abbio.RDS")

summary(fit_abbio)

# Ab: Str affect seedlings (E)
# bf_pol_E   <- bf(multi_pol ~ ConIndex + StrIndex)
# bf_sd_E    <- bf(multi_sd ~ ConIndex + StrIndex)
# bf_seeds_E <- bf(FDSeeds | mi() ~ multi_pol + multi_sd)
# bf_sdlng_E <- bf(FDSdlng | mi() ~ StrIndex)
# 
# fit_ab <- brm(
#   bf_pol_E + bf_sd_E + bf_seeds_E + bf_sdlng_E + set_rescor(FALSE),
#   data = data_scaled, cores = 4, chains = 4, iter = 4000,
#   file = "output/SEM_models/fit_ab"
# )

fit_ab <- readRDS("output/SEM_models/fit_ab.RDS")

summary(fit_ab)

# Null: Nothing affects seedlings (D)
# bf_pol_D   <- bf(multi_pol ~ ConIndex + StrIndex)
# bf_sd_D    <- bf(multi_sd ~ ConIndex + StrIndex)
# bf_seeds_D <- bf(FDSeeds | mi() ~ multi_pol + multi_sd) 
# bf_sdlng_D <- bf(FDSdlng | mi() ~ 1)         
# 
# fit_null <- brm(
#   bf_pol_D + bf_sd_D + bf_seeds_D + bf_sdlng_D + set_rescor(FALSE),
#   data = data_scaled, cores = 4, chains = 4, iter = 4000,
#   file = "output/SEM_models/fit_null"
# )

fit_null <- readRDS("output/SEM_models/fit_null.RDS")

summary(fit_null)

# to allow NAs in model comparison:
# 
# # brms divides the 62 plots in 10 random groups 
# kfold_sat <- kfold(fit_saturated, K = 10, resp = c("multipol", "multisd"), chains = 4, seed = 1212) 
# kfold_med <- kfold(fit_mediated, K = 10, resp = c("multipol", "multisd"), chains = 4, seed = 1212)  
# kfold_abbio <- kfold(fit_abbio, K = 10, resp = c("multipol", "multisd"), chains = 4, seed = 1212) 
# kfold_null <- kfold(fit_null, K = 10, resp = c("multipol", "multisd"), chains = 4, seed = 1212) 

# check if the absolute difference |(elpd_diff)| is two times larger than se (se_diff).
# if not, models are the same
# comp_kfold <- loo_compare(kfold_sat, kfold_med, kfold_abbio, kfold_null)
# comp_kfold

complete_rows <- complete.cases(data_scaled[,c("multi_pol", "multi_sd", "ConIndex", "StrIndex", "FDSeeds", "FDSdlng")]) 
complete_data <- data_scaled[complete_rows,c("multi_pol", "multi_sd", "ConIndex", "StrIndex", "FDSeeds", "FDSdlng")]

loo_sat <- loo(fit_saturated, newdata = complete_data)
loo_med <- loo(fit_mediated, newdata = complete_data)
loo_abbio <- loo(fit_abbio, newdata = complete_data)
loo_ab <- loo(fit_ab, newdata = complete_data)
loo_null <- loo(fit_null, newdata = complete_data)

comp_loo <- loo_compare(loo_sat, loo_med, loo_abbio, loo_ab, loo_null)

## Final diagnostics:

pp_check(fit_ab, resp = "multipol")
pp_check(fit_ab, resp = "multisd")
pp_check(fit_ab, resp = "FDSeeds")
pp_check(fit_ab, resp = "FDSdlng")

bayes_R2(fit_ab, na.rm = T)
conditional_effects(fit_ab, resp = "FDSdlng")

# Função para calcular R2 de uma variável específica em modelos multivariados
get_r2_manual <- function(model, resp_name, sigma_name) {
  # Extrai as predições (linear predictor)
  mu <- fitted(model, resp = resp_name, summary = FALSE)
  var_mu <- apply(mu, 1, var)
  
  # Extrai o sigma do posterior
  sigma <- as_draws_df(model)[[sigma_name]]
  var_epsilon <- sigma^2
  
  # Calcula o R2 para cada draw
  r2_draws <- var_mu / (var_mu + var_epsilon)
  
  # Retorna o resumo (média e intervalos)
  return(data.frame(
    Estimate = mean(r2_draws),
    Q2.5 = quantile(r2_draws, 0.025),
    Q97.5 = quantile(r2_draws, 0.975)
  ))
}

# Exemplos de uso:
r2_seeds <- get_r2_manual(fit_abbio, "FDSeeds", "sigma_FDSeeds")
r2_sdlng <- get_r2_manual(fit_abbio, "FDSdlng", "sigma_FDSdlng")

# checking if it works:
r2_pol <- get_r2_manual(fit_abbio, "multipol", "sigma_multipol")
r2_sd <- get_r2_manual(fit_abbio, "multisd", "sigma_multisd")
