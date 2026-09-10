#### MCMC diagnostics summary ####

library(flextable)
library(officer)

samples <- readRDS("output/model_posteriors_con.rds")
samples_all <- as.mcmc(do.call(rbind, samples))

group_labels <- c(
  "1" = "Bees",
  "2" = "Moths",
  "3" = "Bats (pollination)",
  "4" = "Bats (seed dispersal)",
  "5" = "Birds",
  "6" = "Non-flying mammals"
)

table_long <- bind_rows(
  summarise_param("theta_0",   "A0",   samples, samples_all),
  summarise_param("theta_inf", "AOG",  samples, samples_all),
  summarise_param("alpha_con", "\u03b1",    samples, samples_all),
  summarise_param("beta_con",  "\u03b2",    samples, samples_all),
  summarise_param("sigma_old", "\u03a4OG",  samples, samples_all, transform_fn = tau),
  summarise_param("sigma_rec", "\u03a4REC", samples, samples_all, transform_fn = tau)
)

param_order <- c("A0", "AOG", "\u03b1", "\u03b2", "\u03a4OG", "\u03a4REC")
group_order <- group_labels

table_final <- table_long %>%
  mutate(
    Parameter = factor(Parameter, levels = param_order),
    Group = factor(Group, levels = group_order)
  ) %>%
  arrange(Group, Parameter) %>%
  mutate(
    across(c(Median, `2.5%`, `25%`, `75%`, `97.5%`, PSRF), ~ round(.x, 2)),
    Neff = round(Neff, 1),
    `Parameter ` = paste(Parameter, Group)   # combined label, matching reference table style
  ) %>%
  select(`Parameter ` , Median, `2.5%`, `25%`, `75%`, `97.5%`, Neff, PSRF)

write.csv(table_final, "output/TableS1_diagnostics.csv", row.names = FALSE)
