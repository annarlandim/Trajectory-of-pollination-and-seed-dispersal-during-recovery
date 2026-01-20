## Data

data <- read.csv("data/processed/model_df.csv")

zeros <- data %>% 
  map_lgl(~ any(. == 0, na.rm = TRUE))
print(zeros)

# to avoid problems with log(0)
eps <- 1e-6
data$AbSdlng <- pmax(data$AbSdlng, eps)
data$RichSdlng <- pmax(data$RichSdlng, eps)
data$ShnSdlng <- pmax(data$ShnSdlng, eps)

weights_df <- read.csv("data/processed/weights_df.csv")

par(mfrow = c(1, 3))
plot(FDBees ~ RegTime, data)
plot(FDMoths ~ RegTime, data)
plot(FDBat_pol ~ RegTime, data)

par(mfrow = c(1, 3))
plot(log(FDBats) ~ RegTime, data)
plot(log(FDBirds) ~ RegTime, data)
plot(log(FDNf) ~ RegTime, data)

par(mfrow = c(1, 4))
plot(AbSeeds ~ RegTime, data)
plot(RichSeeds ~ RegTime, data)
plot(ShnSeeds ~ RegTime, data)
plot(FDSeeds ~ RegTime, data)

par(mfrow = c(1, 4))
plot(AbSdlng ~ RegTime, data)
plot(RichSdlng ~ RegTime, data)
plot(ShnSdlng ~ RegTime, data)
plot(FDSdlng ~ RegTime, data)


model_df <- data %>%
  filter(Treatment3 != "old-growth forest") %>%
  select(RegTime, ConIndex, StrIndex, FDBees, FDMoths, FDBat_pol, FDBats, FDBirds, FDNf, FDSeeds, AbSeeds, RichSeeds, ShnSeeds, FDSdlng, AbSdlng, RichSdlng, ShnSdlng) %>%
  mutate(
    logFDBees = log(FDBees),
    logFDMoths = log(FDMoths),
    logFDBat_pol = log(FDBat_pol),
    logFDBats = log(FDBats),
    logFDBirds = log(FDBirds),
    logFDNf = log(FDNf),
    logFDSeeds = log(FDSeeds),
    logAbSeeds = log(AbSeeds),
    logRichSeeds = log(RichSeeds),
    logShnSeeds = log(ShnSeeds),
    logFDSdlng = log(FDSdlng),
    logAbSdlng = log(AbSdlng),
    logRichSdlng = log(RichSdlng),
    logShnSdlng = log(ShnSdlng),
    strc_z = scale(StrIndex)[,1],
    con_z = scale(ConIndex)[,1],
    time_z = scale(RegTime)[,1]
  )

cor.test(model_df$RegTime, model_df$ConIndex)
cor.test(model_df$RegTime, model_df$StrIndex)
cor.test(model_df$StrIndex, model_df$ConIndex)

plot(AbSeeds ~ ConIndex, model_df)
plot(RichSeeds ~ ConIndex, model_df)
plot(ShnSeeds ~ ConIndex, model_df)
plot(FDSeeds ~ ConIndex, model_df)

w_pol <- weights_df[weights_df$index %in% c(1:3), "prctg"] / sum(weights_df[weights_df$index %in% c(1:3), "prctg"])
w_sd <- weights_df[weights_df$index %in% c(4:6), "prctg"] / sum(weights_df[weights_df$index %in% c(4:6), "prctg"])

process_df <- model_df %>%
  rowwise() %>%
  mutate(
    MF_pollination = sum(
      c(FDBees, FDMoths, FDBat_pol) * w_pol,
      na.rm = TRUE
    ),
    MF_dispersal = sum(
      c(FDBats, FDBirds, FDNf) * w_sd,
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  mutate(
    MF_pollination_z = scale(MF_pollination)[,1],
    MF_dispersal_z = scale(MF_dispersal)[,1]
  )

# Univariate models and then, correlation between residuals (maximizes N per group)

## Pollination

cor.test(model_df$logFDBees, model_df$logFDMoths)
cor.test(model_df$logFDBees, model_df$logFDBat_pol)
cor.test(model_df$logFDMoths, model_df$logFDBat_pol)

### Bees

fit_bee_pol_t <- brm(
  logFDBees ~ time_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_bee_pol_t

fit_bee_pol_s <- brm(
  logFDBees ~ strc_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_bee_pol_s

# compare models: the higher the elpd (expected log predictive density), 
# the better the predictive power of the model
# when diff < 2 * SE, models are indistinguishable
loo_bee_t <- loo(fit_bee_pol_t)
loo_bee_s <- loo(fit_bee_pol_s)
loo_compare(loo_bee_t, loo_bee_s)

fit_bee_pol_c <- brm(
  logFDBees ~ con_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_bee_pol_c

### Moths

fit_moth_pol_t <- brm(
  logFDMoths ~ time_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_moth_pol_t

fit_moth_pol_s <- brm(
  logFDMoths ~ strc_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_moth_pol_s

# compare models: the higher the elpd (expected log predictive density), 
# the better the predictive power of the model
# when diff < 2 * SE, models are indistinguishable
loo_moth_t <- loo(fit_moth_pol_t)
loo_moth_s <- loo(fit_moth_pol_s)
loo_compare(loo_moth_t, loo_moth_s)

fit_moth_pol_c <- brm(
  logFDMoths ~ con_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_moth_pol_c

### Bat_pol

fit_bat_pol_t <- brm(
  logFDBat_pol ~ time_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_bat_pol_t

fit_bat_pol_s <- brm(
  logFDBat_pol ~ strc_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_bat_pol_s

# compare models: the higher the elpd (expected log predictive density), 
# the better the predictive power of the model
# when diff < 2 * SE, models are indistinguishable
loo_bat_p_t <- loo(fit_bat_pol_t)
loo_bat_p_s <- loo(fit_bat_pol_s)
loo_compare(loo_bat_p_t, loo_bat_p_s)

fit_bat_pol_c <- brm(
  logFDBat_pol ~ con_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_bat_pol_c

## Seed dispersal

cor.test(model_df$logFDBats, model_df$logFDBirds)
cor.test(model_df$logFDBats, model_df$logFDNf)
cor.test(model_df$logFDBirds, model_df$logFDNf)

### Bees

fit_bat_sd_t <- brm(
  logFDBats ~ time_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_bat_sd_t

fit_bat_sd_s <- brm(
  logFDBats ~ strc_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_bat_sd_s

# compare models: the higher the elpd (expected log predictive density), 
# the better the predictive power of the model
# when |diff| < 2 * SE, models are indistinguishable
loo_bat_sd_t <- loo(fit_bat_sd_t)
loo_bat_sd_s <- loo(fit_bat_sd_s)
loo_compare(loo_bat_sd_t, loo_bat_sd_s)

fit_bat_sd_c <- brm(
  logFDBats ~ con_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_bat_sd_c

### Birds

fit_bird_sd_t <- brm(
  logFDBirds ~ time_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_bird_sd_t

fit_bird_sd_s <- brm(
  logFDBirds ~ strc_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_bird_sd_s

# compare models: the higher the elpd (expected log predictive density), 
# the better the predictive power of the model
# when diff < 2 * SE, models are indistinguishable
loo_bird_t <- loo(fit_bird_sd_t)
loo_bird_s <- loo(fit_bird_sd_s)
loo_compare(loo_bird_t, loo_bird_s)

#k ≤ 0.7 → ok
# 0.7 < k ≤ 1 → problemático
# k > 1 → LOO inválido

fit_bird_sd_c <- brm(
  logFDBirds ~ con_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_bird_sd_c

### Non-flying mammals

fit_nf_sd_t <- brm(
  logFDNf ~ time_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_nf_sd_t

fit_nf_sd_s <- brm(
  logFDNf ~ strc_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_nf_sd_s

# compare models: the higher the elpd (expected log predictive density), 
# the better the predictive power of the model
# when |diff| < 2 * SE, models are indistinguishable
loo_nf_sd_t <- loo(fit_nf_sd_t)
loo_nf_sd_s <- loo(fit_nf_sd_s)
loo_compare(loo_nf_sd_t, loo_nf_sd_s)

fit_nf_sd_c <- brm(
  logFDNf ~ con_z,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_nf_sd_c

### Ploting:

get_slope <- function(model, predictor_name) {
  fe <- fixef(model)
  return(fe[predictor_name, c("Estimate", "Q2.5", "Q97.5")])
}

groups <- list(
  "Bees"    = list(t = fit_bee_pol_t, s = fit_bee_pol_s, c = fit_bee_pol_c, var = "time_z"),
  "Moths"   = list(t = fit_moth_pol_t, s = fit_moth_pol_s, c = fit_moth_pol_c, var = "time_z"),
  "Bat Pol" = list(t = fit_bat_pol_t, s = fit_bat_pol_s, c = fit_bat_pol_c, var = "time_z"),
  "Bats"    = list(t = fit_bat_sd_t,  s = fit_bat_sd_s,  c = fit_bat_sd_c,  var = "time_z"),
  "Birds"   = list(t = fit_bird_sd_t, s = fit_bird_sd_s, c = fit_bird_sd_c, var = "time_z"),
  "NF Mamm" = list(t = fit_nf_sd_t,   s = fit_nf_sd_s,   c = fit_nf_sd_c,   var = "time_z")
)

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

cols <- c("Time" = "#D95F02", "Structure" = "#1B9E77", "Connectivity" = "#7570B3")

for (i in seq_along(groups)) {
  name <- names(groups)[i]
  g <- groups[[i]]
  
  # Extrair valores (Atenção: o nome da variável muda no modelo s e c)
  res_t <- get_slope(g$t, "time_z")
  res_s <- get_slope(g$s, "strc_z")
  res_c <- get_slope(g$c, "con_z")
  
  res_all <- rbind(res_t, res_s, res_c)
  
  # Criar o plot vazio
  plot(NULL, xlim = c(-0.5, 0.5), ylim = c(0.5, 3.5), 
       yaxt = "n", xlab = "Standardized Coefficient (beta)", ylab = "",
       main = name, bty = "n")
  
  # Adicionar linha no zero (sem efeito)
  abline(v = 0, lty = 2, col = "grey50")
  
  # Adicionar os nomes dos preditores no eixo Y
  axis(2, at = 1:3, labels = c("Connectivity", "Structure", "Time"), las = 1)
  
  # Desenhar os intervalos e pontos (de baixo para cima: C, S, T)
  for (j in 1:3) {
    y_pos <- 4 - j # Inverter para Time ficar no topo
    # color <- rev(cols)[j]
    
    # Linha do intervalo de credibilidade
    segments(res_all[j, 2], y_pos, res_all[j, 3], y_pos, col = "black", lwd = 2)
    
    # Ponto da mediana
    # Se o intervalo não cruzar o zero, o ponto é preenchido
    is_sig <- (res_all[j, 2] > 0 | res_all[j, 3] < 0)
    points(res_all[j, 1], y_pos, pch = ifelse(is_sig, 19, 21), 
           bg = "white", col = "black", cex = 1.5)
  }
}

mtext("Drivers of Ecological Recovery by Group", outer = TRUE, side = 3, line = 0, font = 2, cex = 1.2)

# Correlation between residuals:

## Pollination

pol_ctrl_str <- bf(logFDBees ~ strc_z) +
  bf(logFDMoths ~ strc_z) +
  # bf(logFDBat_pol ~ strc_z) +
  set_rescor(TRUE)

fit_pol_ctrl_str <- brm(
  pol_ctrl_str,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)

## Seed dispersal

sd_ctrl_str <- bf(logFDBats ~ strc_z) +
  bf(logFDBirds ~ strc_z) +
  # bf(logFDNf ~ strc_z) +
  set_rescor(TRUE)

fit_sd_ctrl_str <- brm(
  sd_ctrl_str,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)

###################################
########## Effects on plants: #####
###################################

## Seeds ####

fit_seeds_str_fd <- brm(
  logFDSeeds ~ strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_str_fd

fit_seeds_con_fd <- brm(
  logFDSeeds ~ con_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_con_fd

fit_seeds_str_ab <- brm(
  logAbSeeds ~ strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_str_ab

fit_seeds_con_ab <- brm(
  logAbSeeds ~ con_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_con_ab

fit_seeds_str_rich <- brm(
  logRichSeeds ~ strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_str_rich

fit_seeds_con_rich <- brm(
  logRichSeeds ~ con_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_con_rich

fit_seeds_str_shn <- brm(
  logShnSeeds ~ strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_str_shn

fit_seeds_con_shn <- brm(
  logShnSeeds ~ con_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_con_shn

### Pollination:

#### FD

fit_seeds_pol_mf_fd <- brm(
  logFDSeeds ~ MF_pollination_z + strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_pol_mf_fd

fit_seeds_pol_fd <- brm(
  logFDSeeds ~ logFDBees + logFDMoths + strc_z,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_pol_fd

#### Abundance

fit_seeds_pol_mf_ab <- brm(
  logAbSeeds ~ MF_pollination_z + strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_pol_mf_ab

fit_seeds_pol_ab <- brm(
  logAbSeeds ~ logFDBees + logFDMoths + strc_z,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_pol_ab

#### Richness

fit_seeds_pol_mf_rich <- brm(
  logRichSeeds ~ MF_pollination_z + strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_pol_mf_rich

fit_seeds_pol_rich <- brm(
  logRichSeeds ~ logFDBees + logFDMoths + strc_z,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_pol_rich

#### Shannon

fit_seeds_pol_mf_shn <- brm(
  logShnSeeds ~ MF_pollination_z + strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_pol_mf_shn

fit_seeds_pol_shn <- brm(
  logShnSeeds ~ logFDBees + logFDMoths + strc_z,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_pol_shn

### Seed dispersal:

#### FD

fit_seeds_sd_mf_fd <- brm(
  logFDSeeds ~ MF_dispersal_z + strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_sd_mf_fd

fit_seeds_sd_fd <- brm(
  logFDSeeds ~ logFDBats + logFDBirds + strc_z,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_sd_fd

#### Abundance

fit_seeds_sd_mf_ab <- brm(
  logAbSeeds ~ MF_dispersal_z + strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_sd_mf_ab

fit_seeds_sd_ab <- brm(
  logAbSeeds ~ logFDBats + logFDBirds + strc_z,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_sd_ab

#### Richness

fit_seeds_sd_mf_rich <- brm(
  logRichSeeds ~ MF_dispersal_z + strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_sd_mf_rich

fit_seeds_sd_rich <- brm(
  logRichSeeds ~ logFDBats + logFDBirds + strc_z,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_sd_rich

#### Shannon

fit_seeds_sd_mf_shn <- brm(
  logShnSeeds ~ MF_dispersal_z + strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_sd_mf_shn

fit_seeds_sd_shn <- brm(
  logShnSeeds ~ logFDBats + logFDBirds + strc_z,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_seeds_sd_shn

## Seedlings ####

fit_sdlng_str_fd <- brm(
  logFDSdlng ~ strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_str_fd

fit_sdlng_con_fd <- brm(
  logFDSdlng ~ con_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_con_fd

fit_sdlng_str_ab <- brm(
  logAbSdlng ~ strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_str_ab

fit_sdlng_con_ab <- brm(
  logAbSdlng ~ con_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_con_ab

fit_sdlng_str_rich <- brm(
  logRichSdlng ~ strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_str_rich

fit_sdlng_con_rich <- brm(
  logRichSdlng ~ con_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_con_rich

fit_sdlng_str_shn <- brm(
  logShnSdlng ~ strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_str_shn

fit_sdlng_con_shn <- brm(
  logShnSdlng ~ con_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_con_shn

### Pollination:

#### FD

fit_sdlng_pol_mf_fd <- brm(
  logFDSdlng ~ MF_pollination_z + strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_pol_mf_fd

fit_sdlng_pol_fd <- brm(
  logFDSdlng ~ logFDBees + logFDMoths + strc_z,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_pol_fd

#### Abundance

fit_sdlng_pol_mf_ab <- brm(
  logAbSdlng ~ MF_pollination_z + strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_pol_mf_ab

fit_sdlng_pol_ab <- brm(
  logAbSdlng ~ logFDBees + logFDMoths + strc_z,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_pol_ab

#### Richness

fit_sdlng_pol_mf_rich <- brm(
  logRichSdlng ~ MF_pollination_z + strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_pol_mf_rich

fit_sdlng_pol_rich <- brm(
  logRichSdlng ~ logFDBees + logFDMoths + strc_z,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_pol_rich

#### Shannon

fit_sdlng_pol_mf_shn <- brm(
  logShnSdlng ~ MF_pollination_z + strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_pol_mf_shn

fit_sdlng_pol_shn <- brm(
  logShnSdlng ~ logFDBees + logFDMoths + strc_z,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_pol_shn

### Seed dispersal:

#### FD

fit_sdlng_sd_mf_fd <- brm(
  logFDSdlng ~ MF_dispersal_z + strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_sd_mf_fd

fit_sdlng_sd_fd <- brm(
  logFDSdlng ~ logFDBats + logFDBirds + strc_z,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_sd_fd

#### Abundance

fit_sdlng_sd_mf_ab <- brm(
  logAbSdlng ~ MF_dispersal_z + strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_sd_mf_ab

fit_sdlng_sd_ab <- brm(
  logAbSdlng ~ logFDBats + logFDBirds + strc_z,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_sd_ab

#### Richness

fit_sdlng_sd_mf_rich <- brm(
  logRichSdlng ~ MF_dispersal_z + strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_sd_mf_rich

fit_sdlng_sd_rich <- brm(
  logRichSdlng ~ logFDBats + logFDBirds + strc_z,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_sd_rich

#### Shannon

fit_sdlng_sd_mf_shn <- brm(
  logShnSdlng ~ MF_dispersal_z + strc_z,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_sd_mf_shn

fit_sdlng_sd_shn <- brm(
  logShnSdlng ~ logFDBats + logFDBirds + strc_z,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_sd_shn

# Plotting:

## Multifunctionality and env variables

get_beta <- function(model, predictor) {
  if(is.null(model)) return(c(NA, NA, NA))
  f <- fixef(model)
  if(!predictor %in% rownames(f)) return(c(NA, NA, NA))
  return(f[predictor, c("Estimate", "Q2.5", "Q97.5")])
}

metrics <- c("FD", "Abundance", "Richness", "Shannon")
types <- c("Seeds", "Seedlings")

plot_data <- list(
  # Seeds
  "Seeds_FD"   = list(str=fit_seeds_str_fd,   con=fit_seeds_con_fd,   pol=fit_seeds_pol_mf_fd,   sd=fit_seeds_sd_mf_fd),
  "Seeds_Ab"   = list(str=fit_seeds_str_ab,   con=fit_seeds_con_ab,   pol=fit_seeds_pol_mf_ab,   sd=fit_seeds_sd_mf_ab),
  "Seeds_Rich" = list(str=fit_seeds_str_rich, con=fit_seeds_con_rich, pol=fit_seeds_pol_mf_rich, sd=fit_seeds_sd_mf_rich),
  "Seeds_Shn"  = list(str=fit_seeds_str_shn,  con=fit_seeds_con_shn,  pol=fit_seeds_pol_mf_shn,  sd=fit_seeds_sd_mf_shn),
  # Seedlings
  "Sdlng_FD"   = list(str=fit_sdlng_str_fd,   con=fit_sdlng_con_fd,   pol=fit_sdlng_pol_mf_fd,   sd=fit_sdlng_sd_mf_fd),
  "Sdlng_Ab"   = list(str=fit_sdlng_str_ab,   con=fit_sdlng_con_ab,   pol=fit_sdlng_pol_mf_ab,   sd=fit_sdlng_sd_mf_ab),
  "Sdlng_Rich" = list(str=fit_sdlng_str_rich, con=fit_sdlng_con_rich, pol=fit_sdlng_pol_mf_rich, sd=fit_sdlng_sd_mf_rich),
  "Sdlng_Shn"  = list(str=fit_sdlng_str_shn,  con=fit_sdlng_con_shn,  pol=fit_sdlng_pol_mf_shn,  sd=fit_sdlng_sd_mf_shn)
)

# pdf("regeneration_drivers_summary.pdf", width = 12, height = 7)

# Layout: 2 linhas (Sementes e Plântulas) x 4 colunas (Métricas)
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0), mgp = c(2, 0.7, 0))

# Cores consistentes com seus drivers
# cols <- c("Structure" = "#1B9E77", "Connectivity" = "#7570B3",
#           "Pollination (MF)" = "#D95F02", "Dispersal (MF)" = "#E7298A")

for (i in seq_along(plot_data)) {
  name <- names(plot_data)[i]
  models <- plot_data[[i]]
  
  # Extrair Betas
  res_str <- get_beta(models$str, "strc_z")
  res_con <- get_beta(models$con, "con_z")
  res_pol <- get_beta(models$pol, "MF_pollination_z")
  res_sd  <- get_beta(models$sd,  "MF_dispersal_z")
  
  res_all <- rbind(res_sd, res_pol, res_con, res_str)
  rownames(res_all) <- c( "Pollination (MF)", "Seed Dispersal (MF)", "Connectivity", "Structure")
  
  # Criar Plot
  plot(NULL, xlim = c(-0.6, 0.6), ylim = c(0.5, 4.5), 
       yaxt = "n", xlab = "Standardized Coefficient", ylab = "",
       main = gsub("_", " ", name), bty = "n")
  
  abline(v = 0, lty = 2, col = "grey60")
  
  # Apenas o primeiro gráfico de cada linha ganha rótulos no eixo Y
  if (i %% 4 == 1) {
    axis(2, at = 1:4, labels = rownames(res_all), las = 1, cex.axis = 0.9)
  }
  
  for (j in 1:4) {
    if (is.na(res_all[j, 1])) next
    
    # Desenhar Intervalo (95% CI)
    segments(res_all[j, 2], j, res_all[j, 3], j, col = "black", lwd = 2)
    
    # Ponto da Mediana (Preenchido se significativo)
    is_sig <- (res_all[j, 2] > 0 | res_all[j, 3] < 0)
    points(res_all[j, 1], j, pch = 21, bg = ifelse(is_sig, "black", "white"), 
           col = "black", cex = 1.5, lwd = 2)
  }
}

mtext("Drivers of Seed and Seedling Recovery", outer = TRUE, side = 3, line = 0, font = 2, cex = 1.3)
# dev.off()

## Pollination and Seed dispersal groups

biotic_groups <- c("logFDBees", "logFDMoths", "logFDBats", "logFDBirds")
group_labels <- c("Bees", "Moths", "Bats", "Birds")
group_cols <- c("#E41A1C", "#984EA3", "black", "#377EB8") # Vermelho, Roxo, Preto, Azul

biotic_plot_data <- list(
  "Seeds FD"    = list(pol = fit_seeds_pol_fd,  sd = fit_seeds_sd_fd),
  "Seeds Ab"    = list(pol = fit_seeds_pol_ab,  sd = fit_seeds_sd_ab),
  "Seeds Rich"  = list(pol = fit_seeds_pol_rich,sd = fit_seeds_sd_rich),
  "Seeds Shn"   = list(pol = fit_seeds_pol_shn, sd = fit_seeds_sd_shn),
  
  "Sdlng FD"    = list(pol = fit_sdlng_pol_fd, sd = fit_sdlng_sd_fd), # Preencher quando rodar seedlings ~ pollinators
  "Sdlng Ab"    = list(pol = fit_sdlng_pol_ab, sd = fit_sdlng_sd_ab),
  "Sdlng Rich"  = list(pol = fit_sdlng_pol_rich, sd = fit_sdlng_sd_rich),
  "Sdlng Shn"   = list(pol = fit_sdlng_pol_shn, sd = fit_sdlng_sd_shn)
)

# pdf("biotic_responses_full_summary.pdf", width = 14, height = 8)
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1), oma = c(0, 5, 3, 0), mgp = c(2, 0.7, 0))

for (i in seq_along(biotic_plot_data)) {
  name <- names(biotic_plot_data)[i]
  models <- biotic_plot_data[[i]]
  
  # Criar plot vazio
  plot(NULL, xlim = c(-0.8, 0.8), ylim = c(0.5, 4.5), 
       yaxt = "n", xlab = "Coefficient (Beta)", ylab = "",
       main = name, bty = "l", cex.main = 1.2)
  
  abline(v = 0, lty = 2, col = "grey70")
  
  # Rótulos apenas na primeira coluna de cada linha
  if (i %% 4 == 1) {
    axis(2, at = 1:4, labels = rev(group_labels), las = 1, cex.axis = 1.1)
  }
  
  # Se não houver modelos para este painel, pula para o próximo
  if (is.null(models$pol) & is.null(models$sd)) next
  
  # Extrair e plotar cada grupo
  # Ordem no eixo Y: 1=Birds, 2=Bats, 3=Moths, 4=Bees (de baixo para cima)
  
  # Dispersores (extraídos do modelo de dispersão)
  res_bird <- get_beta(models$sd, "logFDBirds")
  res_bat  <- get_beta(models$sd, "logFDBats")
  # Polinizadores (extraídos do modelo de polinização)
  res_moth <- get_beta(models$pol, "logFDMoths")
  res_bee  <- get_beta(models$pol, "logFDBees")
  
  res_matrix <- rbind(res_bird, res_bat, res_moth, res_bee)
  
  for (j in 1:4) {
    if (is.na(res_matrix[j, 1])) next
    
    curr_col <- rev(group_cols)[j]
    
    # Desenhar Intervalos com "bigodes"
    segments(res_matrix[j, 2], j, res_matrix[j, 3], j, col = curr_col, lwd = 2)
    segments(res_matrix[j, 2], j-0.1, res_matrix[j, 2], j+0.1, col = curr_col, lwd = 2)
    segments(res_matrix[j, 3], j-0.1, res_matrix[j, 3], j+0.1, col = curr_col, lwd = 2)
    
    # Ponto da Mediana (Preenchido se 95% CI não cruza o zero)
    is_sig <- (res_matrix[j, 2] > 0 | res_matrix[j, 3] < 0)
    points(res_matrix[j, 1], j, pch = 21, bg = ifelse(is_sig, curr_col, "white"), 
           col = curr_col, cex = 1.8, lwd = 2)
  }
}

mtext("Biotic Drivers of Tropical Forest Regeneration", outer = TRUE, side = 3, line = 0, font = 2, cex = 1.5)
# dev.off()

plot(process_df$MF_dispersal_z ~ process_df$RegTime)
plot(process_df$logAbSeeds ~ process_df$MF_dispersal_z)
plot(process_df$logRichSeeds ~ process_df$MF_dispersal_z)
plot(process_df$logShnSeeds ~ process_df$MF_dispersal_z)



##### Not including 0s for seedlings:



data.2 <- read.csv("data/processed/model_df.csv")

zeros <- data.2 %>% 
  map_lgl(~ any(. == 0, na.rm = TRUE))
print(zeros)

data.2$AbSdlng[data.2$AbSdlng == 0] <- NA
data.2$RichSdlng[data.2$RichSdlng == 0] <- NA
data.2$ShnSdlng[data.2$ShnSdlng == 0] <- NA

par(mfrow = c(1, 4))
plot(AbSdlng ~ RegTime, data.2)
plot(RichSdlng ~ RegTime, data.2)
plot(ShnSdlng ~ RegTime, data.2)
plot(FDSdlng ~ RegTime, data.2)

model_df.2 <- data.2 %>%
  filter(Treatment3 != "old-growth forest") %>%
  select(RegTime, ConIndex, StrIndex, FDBees, FDMoths, FDBat_pol, FDBats, FDBirds, FDNf, FDSeeds, AbSeeds, RichSeeds, ShnSeeds, FDSdlng, AbSdlng, RichSdlng, ShnSdlng) %>%
  mutate(
    logFDBees = log(FDBees),
    logFDMoths = log(FDMoths),
    logFDBat_pol = log(FDBat_pol),
    logFDBats = log(FDBats),
    logFDBirds = log(FDBirds),
    logFDNf = log(FDNf),
    logFDSeeds = log(FDSeeds),
    logAbSeeds = log(AbSeeds),
    logRichSeeds = log(RichSeeds),
    logShnSeeds = log(ShnSeeds),
    logFDSdlng = log(FDSdlng),
    logAbSdlng = log(AbSdlng),
    logRichSdlng = log(RichSdlng),
    logShnSdlng = log(ShnSdlng),
    strc_z = scale(StrIndex)[,1],
    con_z = scale(ConIndex)[,1],
    time_z = scale(RegTime)[,1]
  )
process_df.2 <- model_df.2 %>%
  rowwise() %>%
  mutate(
    MF_pollination = sum(
      c(FDBees, FDMoths, FDBat_pol) * w_pol,
      na.rm = TRUE
    ),
    MF_dispersal = sum(
      c(FDBats, FDBirds, FDNf) * w_sd,
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  mutate(
    MF_pollination_z = scale(MF_pollination)[,1],
    MF_dispersal_z = scale(MF_dispersal)[,1]
  )

fit_sdlng_str_fd.2 <- brm(
  logFDSdlng ~ strc_z,
  data = process_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_str_fd.2

fit_sdlng_con_fd.2 <- brm(
  logFDSdlng ~ con_z,
  data = process_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_con_fd.2

fit_sdlng_str_ab.2 <- brm(
  logAbSdlng ~ strc_z,
  data = process_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_str_ab.2

fit_sdlng_con_ab.2 <- brm(
  logAbSdlng ~ con_z,
  data = process_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_con_ab.2

fit_sdlng_str_rich.2 <- brm(
  logRichSdlng ~ strc_z,
  data = process_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_str_rich.2

fit_sdlng_con_rich.2 <- brm(
  logRichSdlng ~ con_z,
  data = process_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_con_rich.2

fit_sdlng_str_shn.2 <- brm(
  logShnSdlng ~ strc_z,
  data = process_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_str_shn.2

fit_sdlng_con_shn.2 <- brm(
  logShnSdlng ~ con_z,
  data = process_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_con_shn.2

### Pollination:

#### FD

fit_sdlng_pol_mf_fd.2 <- brm(
  logFDSdlng ~ MF_pollination_z + strc_z,
  data = process_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_pol_mf_fd.2

fit_sdlng_pol_fd.2 <- brm(
  logFDSdlng ~ logFDBees + logFDMoths + strc_z,
  data = model_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_pol_fd.2

#### Abundance

fit_sdlng_pol_mf_ab.2 <- brm(
  logAbSdlng ~ MF_pollination_z + strc_z,
  data = process_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_pol_mf_ab.2

fit_sdlng_pol_ab.2 <- brm(
  logAbSdlng ~ logFDBees + logFDMoths + strc_z,
  data = model_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_pol_ab.2

#### Richness

fit_sdlng_pol_mf_rich.2 <- brm(
  logRichSdlng ~ MF_pollination_z + strc_z,
  data = process_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_pol_mf_rich.2

fit_sdlng_pol_rich.2 <- brm(
  logRichSdlng ~ logFDBees + logFDMoths + strc_z,
  data = model_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_pol_rich.2

#### Shannon

fit_sdlng_pol_mf_shn.2 <- brm(
  logShnSdlng ~ MF_pollination_z + strc_z,
  data = process_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_pol_mf_shn.2

fit_sdlng_pol_shn.2 <- brm(
  logShnSdlng ~ logFDBees + logFDMoths + strc_z,
  data = model_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_pol_shn.2

### Seed dispersal:

#### FD

fit_sdlng_sd_mf_fd.2 <- brm(
  logFDSdlng ~ MF_dispersal_z + strc_z,
  data = process_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_sd_mf_fd.2

fit_sdlng_sd_fd.2 <- brm(
  logFDSdlng ~ logFDBats + logFDBirds + strc_z,
  data = model_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_sd_fd.2

#### Abundance

fit_sdlng_sd_mf_ab.2 <- brm(
  logAbSdlng ~ MF_dispersal_z + strc_z,
  data = process_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_sd_mf_ab.2

fit_sdlng_sd_ab.2 <- brm(
  logAbSdlng ~ logFDBats + logFDBirds + strc_z,
  data = model_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_sd_ab.2

#### Richness

fit_sdlng_sd_mf_rich.2 <- brm(
  logRichSdlng ~ MF_dispersal_z + strc_z,
  data = process_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_sd_mf_rich.2

fit_sdlng_sd_rich.2 <- brm(
  logRichSdlng ~ logFDBats + logFDBirds + strc_z,
  data = model_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_sd_rich.2

#### Shannon

fit_sdlng_sd_mf_shn.2 <- brm(
  logShnSdlng ~ MF_dispersal_z + strc_z,
  data = process_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_sd_mf_shn.2

fit_sdlng_sd_shn.2 <- brm(
  logShnSdlng ~ logFDBats + logFDBirds + strc_z,
  data = model_df.2,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sdlng_sd_shn.2

## Multifunctionality and env variables

get_beta <- function(model, predictor) {
  if(is.null(model)) return(c(NA, NA, NA))
  f <- fixef(model)
  if(!predictor %in% rownames(f)) return(c(NA, NA, NA))
  return(f[predictor, c("Estimate", "Q2.5", "Q97.5")])
}

metrics <- c("FD", "Abundance", "Richness", "Shannon")
types <- c("Seedling")

plot_data <- list(
  # Seedlings
  "Sdlng_FD"   = list(str=fit_sdlng_str_fd,   con=fit_sdlng_con_fd,   pol=fit_sdlng_pol_mf_fd,   sd=fit_sdlng_sd_mf_fd),
  "Sdlng_Ab"   = list(str=fit_sdlng_str_ab.2,   con=fit_sdlng_con_ab.2,   pol=fit_sdlng_pol_mf_ab.2,   sd=fit_sdlng_sd_mf_ab.2),
  "Sdlng_Rich" = list(str=fit_sdlng_str_rich.2, con=fit_sdlng_con_rich.2, pol=fit_sdlng_pol_mf_rich.2, sd=fit_sdlng_sd_mf_rich.2),
  "Sdlng_Shn"  = list(str=fit_sdlng_str_shn.2,  con=fit_sdlng_con_shn.2,  pol=fit_sdlng_pol_mf_shn.2,  sd=fit_sdlng_sd_mf_shn.2)
)

# pdf("regeneration_drivers_summary.pdf", width = 12, height = 7)

# Layout: 2 linhas (Sementes e Plântulas) x 4 colunas (Métricas)
par(mfrow = c(1, 4), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0), mgp = c(2, 0.7, 0))

# Cores consistentes com seus drivers
# cols <- c("Structure" = "#1B9E77", "Connectivity" = "#7570B3",
#           "Pollination (MF)" = "#D95F02", "Dispersal (MF)" = "#E7298A")

for (i in seq_along(plot_data)) {
  name <- names(plot_data)[i]
  models <- plot_data[[i]]
  
  # Extrair Betas
  res_str <- get_beta(models$str, "strc_z")
  res_con <- get_beta(models$con, "con_z")
  res_pol <- get_beta(models$pol, "MF_pollination_z")
  res_sd  <- get_beta(models$sd,  "MF_dispersal_z")
  
  res_all <- rbind(res_sd, res_pol, res_con, res_str)
  rownames(res_all) <- c( "Pollination (MF)", "Seed Dispersal (MF)", "Connectivity", "Structure")
  
  # Criar Plot
  plot(NULL, xlim = c(-1, 1), ylim = c(0.5, 4.5), 
       yaxt = "n", xlab = "Standardized Coefficient", ylab = "",
       main = gsub("_", " ", name), bty = "n")
  
  abline(v = 0, lty = 2, col = "grey60")
  
  # Apenas o primeiro gráfico de cada linha ganha rótulos no eixo Y
  if (i %% 4 == 1) {
    axis(2, at = 1:4, labels = rownames(res_all), las = 1, cex.axis = 0.9)
  }
  
  for (j in 1:4) {
    if (is.na(res_all[j, 1])) next
    
    # Desenhar Intervalo (95% CI)
    segments(res_all[j, 2], j, res_all[j, 3], j, col = "black", lwd = 2)
    
    # Ponto da Mediana (Preenchido se significativo)
    is_sig <- (res_all[j, 2] > 0 | res_all[j, 3] < 0)
    points(res_all[j, 1], j, pch = 21, bg = ifelse(is_sig, "black", "white"), 
           col = "black", cex = 1.5, lwd = 2)
  }
}

mtext("Drivers of Seed and Seedling Recovery", outer = TRUE, side = 3, line = 0, font = 2, cex = 1.3)
# dev.off()

## Pollination and Seed dispersal groups

biotic_groups <- c("logFDBees", "logFDMoths", "logFDBats", "logFDBirds")
group_labels <- c("Bees", "Moths", "Bats", "Birds")
group_cols <- c("#E41A1C", "#984EA3", "black", "#377EB8") # Vermelho, Roxo, Preto, Azul

biotic_plot_data <- list(
  "Sdlng FD"    = list(pol = fit_sdlng_pol_fd, sd = fit_sdlng_sd_fd), # Preencher quando rodar seedlings ~ pollinators
  "Sdlng Ab"    = list(pol = fit_sdlng_pol_ab.2, sd = fit_sdlng_sd_ab.2),
  "Sdlng Rich"  = list(pol = fit_sdlng_pol_rich.2, sd = fit_sdlng_sd_rich.2),
  "Sdlng Shn"   = list(pol = fit_sdlng_pol_shn.2, sd = fit_sdlng_sd_shn.2)
)

# pdf("biotic_responses_full_summary.pdf", width = 14, height = 8)
par(mfrow = c(1, 4), mar = c(4, 4, 3, 1), oma = c(0, 5, 3, 0), mgp = c(2, 0.7, 0))

for (i in seq_along(biotic_plot_data)) {
  name <- names(biotic_plot_data)[i]
  models <- biotic_plot_data[[i]]
  
  # Criar plot vazio
  plot(NULL, xlim = c(-1.5, 1.5), ylim = c(0.5, 4.5), 
       yaxt = "n", xlab = "Coefficient (Beta)", ylab = "",
       main = name, bty = "l", cex.main = 1.2)
  
  abline(v = 0, lty = 2, col = "grey70")
  
  # Rótulos apenas na primeira coluna de cada linha
  if (i %% 4 == 1) {
    axis(2, at = 1:4, labels = rev(group_labels), las = 1, cex.axis = 1.1)
  }
  
  # Se não houver modelos para este painel, pula para o próximo
  if (is.null(models$pol) & is.null(models$sd)) next
  
  # Extrair e plotar cada grupo
  # Ordem no eixo Y: 1=Birds, 2=Bats, 3=Moths, 4=Bees (de baixo para cima)
  
  # Dispersores (extraídos do modelo de dispersão)
  res_bird <- get_beta(models$sd, "logFDBirds")
  res_bat  <- get_beta(models$sd, "logFDBats")
  # Polinizadores (extraídos do modelo de polinização)
  res_moth <- get_beta(models$pol, "logFDMoths")
  res_bee  <- get_beta(models$pol, "logFDBees")
  
  res_matrix <- rbind(res_bird, res_bat, res_moth, res_bee)
  
  for (j in 1:4) {
    if (is.na(res_matrix[j, 1])) next
    
    curr_col <- rev(group_cols)[j]
    
    # Desenhar Intervalos com "bigodes"
    segments(res_matrix[j, 2], j, res_matrix[j, 3], j, col = curr_col, lwd = 2)
    segments(res_matrix[j, 2], j-0.1, res_matrix[j, 2], j+0.1, col = curr_col, lwd = 2)
    segments(res_matrix[j, 3], j-0.1, res_matrix[j, 3], j+0.1, col = curr_col, lwd = 2)
    
    # Ponto da Mediana (Preenchido se 95% CI não cruza o zero)
    is_sig <- (res_matrix[j, 2] > 0 | res_matrix[j, 3] < 0)
    points(res_matrix[j, 1], j, pch = 21, bg = ifelse(is_sig, curr_col, "white"), 
           col = curr_col, cex = 1.8, lwd = 2)
  }
}

mtext("Biotic Drivers of Tropical Forest Regeneration", outer = TRUE, side = 3, line = 0, font = 2, cex = 1.5)
# dev.off()
###################################################################################
## Multivariate model

# because we are modelling more than one response variable at the same time, allowing for statistical dependency between them

# For each plot i and group j:
# Y_ij ~ LogNormal(mu_ij, sigma_j)
# mu_ij = f_j(RegTime_i) + B_jS * Structure_i + B_jC * Connectivity_i
# ε_i = (ε_i1, ε_i2,...) ~ N(0, Σ), where ε_i is everything that the model does not explain for a certain variable in plot i

# In the model accounting only for recovery time, one assumes that errors of each group in a plot are independent of each other.
# When we say that ε_i ~ N(0, Σ), we say that errors of different groups in the same plot come from a joint distribution. 
# This means that each plot generates a vetor of errors (ε_i1, ε_i2,...). This vector has mean 0 and a covariance matrix Σ 
# Σ is a correlation matrix between the errors of all variables. its diagonal shows the residuals variance and outside
# the diagnoal are the residuals covariances. 
# saying that variables share their error structure = saying their residuals covary after controlling for predictors

# when rescor = FALSE we say that the residuals are independent
# when rescor = TRUE we say that the residuals are correlated



### Using all groups
## 1. Comparing within processes:

# The models below allow me to see if the variables modelled share their errors structure
# set_rescor(TRUE) checks if the residuals of the models are correlated

# Pollination

cor.test(model_df$logFDBees, model_df$logFDMoths)

pol_t <- bf(logFDBees ~ time_z) +
  bf(logFDMoths ~ time_z) +
  # bf(logFDBat_pol ~ time_z) +
  set_rescor(TRUE)

fit_pol_t <- brm(
  pol_t,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_pol_t
# How to interpret:
# Intercept: if the interval does not cross 0, the estimate is good
# in log-scale an intercept is exp(intercept), so, for instance exp(0.23) = 1.26, meaning variable is 26% above reference value
# Slope: if -0,08[-0.14, -0.02], this means that for each +1 SD in time, bees' FD decreases ~8% in logscale
# since the CI does not crosses 0, this effect is supported
# If -0.08[-0.2, 0.04], effect is similar, but there is a great uncertainty, so no clear evidence of effect

# sigma give the residuals sd, after removing the effect of time
# so it shows how much the data varies after removing the effect of time
# moths (0.39) have much more variability than bees (0.20)

# residual correlations: when 0.16[-0.16, 0.45] large interval, that includes 0, no strong evidation of residual dependency

pol_s <- bf(logFDBees ~ strc_z) +
  bf(logFDMoths ~ strc_z) +
  # bf(logFDBat_pol ~ strc_z) +
  set_rescor(TRUE)

fit_pol_s <- brm(
  pol_s,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_pol_s

pol_c <- bf(logFDBees ~ con_z) +
  bf(logFDMoths ~ con_z) +
  # bf(logFDBat_pol ~ con_z) +
  set_rescor(TRUE)

fit_pol_c <- brm(
  pol_c,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_pol_c

loo(fit_pol_s, fit_pol_t)

# doesnt really make sense since only structure and time seem to affect pollinators (and these are correlated)
pol_all <- bf(logFDBees ~ time_z + strc_z + con_z) +
  bf(FDMoths ~ time_z + strc_z + con_z) +
  bf(FDBat_pol ~ time_z + strc_z + con_z) +
  set_rescor(TRUE)

fit_pol_all <- brm(
  pol_all,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)

pol_ctrl_str <- bf(logFDBees ~ time_z + strc_z) +
  bf(FDMoths ~ time_z + strc_z) +
  bf(FDBat_pol ~ time_z + strc_z) +
  set_rescor(TRUE)

fit_pol_ctrl_str <- brm(
  pol_ctrl_str,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)

pol_ctrl_con <- bf(logFDBees ~ time_z + con_z) +
  bf(FDMoths ~ time_z + con_z) +
  bf(FDBat_pol ~ time_z + con_z) +
  set_rescor(TRUE)

fit_pol_ctrl_con <- brm(
  pol_ctrl_con,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)

# Seed-dispersal
sd_t <- bf(logFDBats ~ time_z) +
  bf(logFDBirds ~ time_z) +
  # bf(logFDNf ~ time_z) +
  set_rescor(TRUE)

fit_sd_t <- brm(
  sd_t,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_sd_t

sd_s <- bf(logFDBats ~ strc_z) +
  bf(logFDBirds ~ strc_z) +
  # bf(logFDNf ~ time_z) +
  set_rescor(TRUE)

fit_sd_s <- brm(
  sd_s,
  data = model_df,
  family = gaussian(), # I would use the lognormal to be equivalent to the recovery model. However, stan only allows gaussian
  chains = 4, cores = 4, iter = 4000
)
fit_sd_s

sd_sc <- bf(logFDBats ~  strc_z + con_z) +
  bf(logFDBirds ~  strc_z + con_z) +
  # bf(logFDNf ~  strc_z + con_z) +
  set_rescor(TRUE)

fit_sd_sc <- brm(
  sd_sc,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
fit_sd_sc

sd_ctrl_all <- bf(logFDBats ~ time_z + strc_z + con_z) +
  bf(logFDBirds ~ time_z + strc_z + con_z) +
  bf(logFDNf ~ time_z + strc_z + con_z) +
  set_rescor(TRUE)

fit_sd_ctrl_all <- brm(
  sd_ctrl_all,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)

sd_ctrl_str <- bf(logFDBats ~ time_z + strc_z) +
  bf(logFDBirds ~ time_z + strc_z) +
  bf(logFDNf ~ time_z + strc_z) +
  set_rescor(TRUE)

fit_sd_ctrl_str <- brm(
  sd_ctrl_str,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)

sd_ctrl_con <- bf(logFDBats ~ time_z + con_z) +
  bf(logFDBirds ~ time_z + con_z) +
  bf(logFDNf ~ time_z + con_z) +
  set_rescor(TRUE)

fit_sd_ctrl_con <- brm(
  sd_ctrl_con,
  data = model_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)

# comparisons
posterior_summary(fit_sd_raw, pars = "^rescor")
posterior_summary(fit_sd_ctrl_all, pars = "^rescor")
posterior_summary(fit_sd_ctrl_str, pars = "^rescor")
posterior_summary(fit_sd_ctrl_con, pars = "^rescor")

## Effect of pollination and seed dispersal in seeds and seedlings

w_pol <- weights_df[weights_df$index %in% c(1:3), "prctg"] / sum(weights_df[weights_df$index %in% c(1:3), "prctg"])
w_sd <- weights_df[weights_df$index %in% c(4:6), "prctg"] / sum(weights_df[weights_df$index %in% c(4:6), "prctg"])

process_df <- model_df %>%
  rowwise() %>%
  mutate(
    MF_pollination = sum(
      c(FDBees, FDMoths, FDBat_pol) * w_pol,
      na.rm = TRUE
    ),
    MF_dispersal = sum(
      c(FDBats, FDBirds, FDNf) * w_sd,
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  mutate(
    MF_pollination_z = scale(MF_pollination),
    MF_dispersal_z = scale(MF_dispersal)
  )

m_process <- 
  bf(logFDSeeds ~ MF_pollination_z + MF_dispersal_z +
       time_z + strc_z + con_z) +
  bf(logFDSdlng ~ MF_pollination_z + MF_dispersal_z +
       time_z + strc_z + con_z) +
  set_rescor(FALSE)

fit_process <- brm(
  m_process,
  data = process_df,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000
)
posterior_summary(fit_plants, pars = "^rescor")
