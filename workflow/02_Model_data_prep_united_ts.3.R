
#### Data frames ####

master <- read.csv(file = "data/raw/Ecuador Plots_Master_V2.csv")

# Pollination:
int_bees <- read.csv(file = "data/processed/int_bees.csv")
int_moths <- read.csv(file = "data/processed/int_moths.csv")
int_bat_pol <- read.csv(file = "data/processed/int_bat_pol.csv")

# Seed dispersal
int_bats <- read.csv(file = "data/processed/int_bats.csv")
int_birds <- read.csv(file = "data/processed/int_birds.csv")
int_nf <- read.csv(file = "data/processed/int_nf.csv")

# Seeds
seeds <- read.csv(file = "data/processed/seeds.csv")

# Seedlings
seedlings <- read.csv(file = "data/processed/seedlings_ind.csv")

# Traits

traits_bees <- read.csv(file = "data/processed/traits_bees.csv")
traits_moths <- read.csv(file = "data/processed/traits_moth.csv")

traits_bats <- read.csv(file = "data/processed/traits_bats.csv")
traits_birds <- read.csv(file = "data/processed/traits_birds.csv")
traits_nf <- read.csv(file = "data/processed/traits_nf.csv")

traits_plants <- read.csv(file = "data/processed/traits_plants.csv")
# for some reason Eschweilera_rimbachii is duplicated, so:
traits_plants <- traits_plants %>% distinct(plant_species, .keep_all = TRUE)
traits_seeds <- read.csv(file = "data/processed/traits_seeds.csv")
traits_seedlings <- read.csv(file = "data/processed/traits_seedlings.csv")

#### Plants per group:

# Pollination

plants_bees <- unique(int_bees$plant_species)
plants_moths <- unique(int_moths$plant_species)
plants_bats_pol <- unique(int_bat_pol$plant_species)

plants_pol_all <- unique(c(plants_bees, plants_moths, plants_bats_pol))

p_plants_bees <- length(plants_bees)/length(plants_pol_all)
p_plants_moths <- length(plants_moths)/length(plants_pol_all)
p_plants_bats_pol <- length(plants_bats_pol)/length(plants_pol_all)

weights_df_pol <- data.frame(index = 1:3,
                          group = c("Bees", "Moths", "Bats_pol"), 
                          prctg = c(p_plants_bees, p_plants_moths, p_plants_bats_pol))

# Seed dispersal

plants_bats <- unique(int_bats$plant_species)
plants_birds <- unique(int_birds$plant_species)
plants_nf <- unique(int_nf$plant_species)

plants_sd_all <- unique(c(plants_bats, plants_birds, plants_nf))

p_plants_bats <- length(plants_bats)/length(plants_sd_all)
p_plants_birds <- length(plants_birds)/length(plants_sd_all)
p_plants_nf <- length(plants_nf)/length(plants_sd_all)

weights_df_sd <- data.frame( index = 4:6,
                          group = c("Bats", "Birds", "NF"), 
                          prctg = c(p_plants_bats, p_plants_birds, p_plants_nf))

weights_df <- bind_rows(weights_df_pol, weights_df_sd)

# write.csv(weights_df, file = here::here("data", "processed", "weights_df.csv"), row.names = F)

#### Pollination 

int_bees_hbt <- int_bees %>%
  mutate(interaction = paste(plant_species, animal_species, sep = ".")) %>%
  group_by(interaction, Treatment3) %>%
  summarise(Treatment3 = first(Treatment3), animal_species = first(animal_species),
            plant_species = first(plant_species), interaction = first(interaction), 
            .groups = 'drop') %>%
  left_join(traits_bees, by = "animal_species", relationship = "many-to-many") %>%
  left_join(traits_plants, by = "plant_species", relationship = "many-to-many") %>%
  mutate(LogPL = log(prLength), LogWL = log(wingLength), LogH = log(Height), LogCR = log(CorLength), group = "Bees") %>%
  mutate(    
    StdPL = as.numeric(scale(LogPL)),
    StdWL = as.numeric(scale(LogWL)),
    StdH = as.numeric(scale(LogH)),
    StdCR = as.numeric(scale(LogCR))) %>%
  select(Treatment3, interaction, LogPL, LogWL, LogH, LogCR, StdPL, StdWL, StdH, StdCR, group)

int_moths_hbt <- int_moths %>%
  mutate(interaction = paste(plant_species, animal_species, sep = ".")) %>%
  group_by(interaction, Treatment3) %>%
  summarise(Treatment3 = first(Treatment3), animal_species = first(animal_species),
            plant_species = first(plant_species), interaction = first(interaction), 
            .groups = 'drop') %>%
  left_join(traits_moths, by = "animal_species", relationship = "many-to-many") %>%
  left_join(traits_plants, by = "plant_species", relationship = "many-to-many") %>%
  mutate(LogPL = log(prLength), LogWL = log(wingLength), LogH = log(Height), LogCR = log(CorLength), group = "Moths") %>%
  mutate(    
    StdPL = as.numeric(scale(LogPL)),
    StdWL = as.numeric(scale(LogWL)),
    StdH = as.numeric(scale(LogH)),
    StdCR = as.numeric(scale(LogCR))) %>%
  select(Treatment3, interaction, LogPL, LogWL, LogH, LogCR, StdPL, StdWL, StdH, StdCR, group)

int_bat_pol_hbt <- int_bat_pol %>%
  mutate(interaction = paste(plant_species, animal_species, sep = ".")) %>%
  group_by(interaction, Treatment3) %>%
  summarise(Treatment3 = first(Treatment3), animal_species = first(animal_species),
            plant_species = first(plant_species), interaction = first(interaction), 
            .groups = 'drop') %>%
  left_join(traits_bats, by = "animal_species", relationship = "many-to-many") %>%
  left_join(traits_plants, by = "plant_species", relationship = "many-to-many") %>%
  mutate(LogPL = log(GapeWidth), LogH = log(Height), LogCR = log(CorLength), group = "Bat_pol") %>%
  mutate(    
    StdPL = as.numeric(scale(LogPL)),
    StdWL = as.numeric(scale(HWIndex)),
    StdH = as.numeric(scale(LogH)),
    StdCR = as.numeric(scale(LogCR))) %>%
  select(Treatment3, interaction, LogPL, LogWL = HWIndex, LogH, LogCR, StdPL, StdWL, StdH, StdCR, group)

int_bees_plot <- int_bees %>%
  mutate(interaction = paste(plant_species, animal_species, sep = ".")) %>%
  group_by(interaction, Treatment3) %>%
  summarise(Treatment3 = first(Treatment3),
            Plot_ID = first(Plot_ID), RegTime = first(RegTime),
            interaction = first(interaction),
            .groups = "drop") 

int_moths_plot <- int_moths %>%
  mutate(interaction = paste(plant_species, animal_species, sep = ".")) %>%
  group_by(interaction, Treatment3) %>%
  summarise(Treatment3 = first(Treatment3),
            Plot_ID = first(Plot_ID), RegTime = first(RegTime),
            interaction = first(interaction),
            .groups = "drop") 

int_bat_pol_plot <- int_bat_pol %>%
  mutate(interaction = paste(plant_species, animal_species, sep = ".")) %>%
  group_by(interaction, Treatment3) %>%
  summarise(Treatment3 = first(Treatment3),
            Plot_ID = first(Plot_ID), RegTime = first(RegTime),
            interaction = first(interaction),
            .groups = "drop") 

int_pol_hbt <- bind_rows(int_bees_hbt, int_moths_hbt, int_bat_pol_hbt) %>%
  mutate(    
    StdPL = as.numeric(scale(StdPL)),
    StdWL = as.numeric(scale(StdWL)),
    StdH = as.numeric(scale(StdH)),
    StdCR = as.numeric(scale(StdCR))) %>%
  select(Treatment3, interaction, LogPL, LogWL, LogH, LogCR, StdPL, StdWL, StdH, StdCR, group)


#### Seed-dispersal 

int_bats_hbt <- int_bats %>%
  mutate(interaction = paste(plant_species, animal_species, sep = ".")) %>%
  group_by(interaction, Treatment3) %>%
  filter(!plant_species %in% c("Theobroma_cacao", "Manihot_esculenta", "Artocarpus_heterophyllus", "Psidium_guajava",
                               "Borojoa_sp.", "Persea_americana")) %>%
  summarise(Treatment3 = first(Treatment3), animal_species = first(animal_species),
            plant_species = first(plant_species), interaction = first(interaction), 
            .groups = 'drop') %>%
  left_join(traits_bats, by = "animal_species", relationship = "many-to-many") %>%
  left_join(traits_plants, by = "plant_species", relationship = "many-to-many") %>%
  mutate(LogBM = log(BodyMass), LogCM = log(CropMass), LogGW = log(GapeWidth),LogFW = log(FruitWidth), LogH = log(Height), group = "Bats") %>%
  mutate(    
    StdBM = as.numeric(scale(LogBM)),
    StdCM = as.numeric(scale(LogCM)),
    StdGW = as.numeric(scale(LogGW)),
    StdFW = as.numeric(scale(LogFW)),
    StdHWI = as.numeric(scale(HWIndex)),
    StdH = as.numeric(scale(LogH))) %>%
  select(Treatment3, interaction, LogBM, LogCM, LogGW, LogFW, HWIndex, LogH, StdBM, StdCM, StdGW, StdFW, StdHWI, StdH, group)

int_birds_hbt <- int_birds %>%
  mutate(interaction = paste(plant_species, animal_species, sep = ".")) %>%
  group_by(interaction, Treatment3) %>%
  filter(!plant_species %in% c("Theobroma_cacao", "Manihot_esculenta", "Artocarpus_heterophyllus", "Psidium_guajava",
                               "Borojoa_sp.", "Persea_americana")) %>%
  summarise(Treatment3 = first(Treatment3), animal_species = first(animal_species),
            plant_species = first(plant_species), interaction = first(interaction), 
            .groups = 'drop') %>%
  left_join(traits_birds, by = "animal_species", relationship = "many-to-many") %>%
  left_join(traits_plants, by = "plant_species", relationship = "many-to-many") %>%
  mutate(LogBM = log(BodyMass), LogCM = log(CropMass), LogGW = log(BeakWidth), LogFW = log(FruitWidth), LogH = log(Height), group = "Birds") %>%
  mutate(    
    StdBM = as.numeric(scale(LogBM)),
    StdCM = as.numeric(scale(LogCM)),
    StdGW = as.numeric(scale(LogGW)),
    StdFW = as.numeric(scale(LogFW)),
    StdHWI = as.numeric(scale(HWIndex)),
    StdH = as.numeric(scale(LogH))) %>%
  select(Treatment3, interaction, LogBM, LogCM, LogGW, LogFW, HWIndex, LogH, StdBM, StdCM, StdGW, StdFW, StdHWI, StdH, group)

int_nf_hbt <- int_nf %>%
  mutate(interaction = paste(plant_species, animal_species, sep = ".")) %>%
  group_by(interaction, Treatment3) %>%
  filter(!plant_species %in% c("Theobroma_cacao", "Manihot_esculenta", "Artocarpus_heterophyllus", "Psidium_guajava",
                               "Borojoa_sp.", "Persea_americana")) %>%
  summarise(Treatment3 = first(Treatment3), animal_species = first(animal_species),
            plant_species = first(plant_species), interaction = first(interaction), 
            .groups = 'drop') %>%
  left_join(traits_nf, by = "animal_species", relationship = "many-to-many") %>%
  left_join(traits_plants, by = "plant_species", relationship = "many-to-many") %>%
  mutate(LogBM = log(BodyMass), LogCM = log(CropMass), LogGW = log(GapeWidth), LogFW = log(FruitWidth), LogH = log(Height), group = "NF") %>%
  mutate(    
    StdBM = as.numeric(scale(LogBM)),
    StdCM = as.numeric(scale(LogCM)),
    StdGW = as.numeric(scale(LogGW)),
    StdFW = as.numeric(scale(LogFW)),
    StdHWI = NA,
    StdH = as.numeric(scale(LogH))) %>%
  select(Treatment3, interaction, LogBM, LogCM, LogGW, LogFW, LogH, StdBM, StdCM, StdGW, StdFW, StdHWI, StdH, group)

int_bats_plot <- int_bats %>%
  mutate(interaction = paste(plant_species, animal_species, sep = ".")) %>%
  group_by(interaction, Treatment3) %>%
  filter(!plant_species %in% c("Theobroma_cacao", "Manihot_esculenta", "Artocarpus_heterophyllus", "Psidium_guajava",
                               "Borojoa_sp.", "Persea_americana")) %>%
  summarise(Treatment3 = first(Treatment3),
            Plot_ID = first(Plot_ID), RegTime = first(RegTime),
            interaction = first(interaction),
            .groups = "drop") 

int_birds_plot <- int_birds %>%
  mutate(interaction = paste(plant_species, animal_species, sep = ".")) %>%
  group_by(interaction, Treatment3) %>%
  filter(!plant_species %in% c("Theobroma_cacao", "Manihot_esculenta", "Artocarpus_heterophyllus", "Psidium_guajava",
                               "Borojoa_sp.", "Persea_americana")) %>%
  summarise(Treatment3 = first(Treatment3),
            Plot_ID = first(Plot_ID), RegTime = first(RegTime),
            interaction = first(interaction),
            .groups = "drop") 

int_nf_plot <- int_nf %>%
  mutate(interaction = paste(plant_species, animal_species, sep = ".")) %>%
  group_by(interaction, Treatment3) %>%
  filter(!plant_species %in% c("Theobroma_cacao", "Manihot_esculenta", "Artocarpus_heterophyllus", "Psidium_guajava",
                               "Borojoa_sp.", "Persea_americana")) %>%
  summarise(Treatment3 = first(Treatment3),
            Plot_ID = first(Plot_ID), RegTime = first(RegTime),
            interaction = first(interaction),
            .groups = "drop") 

int_sd_hbt <- bind_rows(int_bats_hbt, int_birds_hbt, int_nf_hbt) %>%
  mutate(    
    StdBM = as.numeric(scale(StdBM)),
    StdCM = as.numeric(scale(StdCM)),
    StdGW = as.numeric(scale(StdGW)),
    StdFW = as.numeric(scale(StdFW)),
    StdHWI = as.numeric(scale(StdHWI)),
    StdH = as.numeric(scale(StdH))) %>%
  select(Treatment3, interaction, LogBM, LogCM, LogGW, LogFW, HWIndex, LogH, StdBM, StdCM, StdGW, StdFW, StdHWI, StdH, group)


#### Seeds and seedlings ####

### Seeds

seeds_hbt <- seeds %>%
  group_by(plant_species, Treatment3) %>%
  filter(!plant_species %in% c("Theobroma_cacao", "Manihot_esculenta", "Artocarpus_heterophyllus", "Psidium_guajava",
                              "Borojoa_sp.", "Persea_americana")) %>%
  summarise(Treatment3 = first(Treatment3), 
            Species_name = first(plant_species), 
            .groups = 'drop') %>%
  left_join(traits_seeds, by = "plant_species", relationship = "many-to-many") %>%
  ungroup()%>%
  mutate(LogWidth = log(SeedWidth), LogLength = log(SeedLength), LogWeight = log(SeedWeight)) %>%
  mutate(    
    StdWidth = as.numeric(scale(LogWidth)),
    StdLength = as.numeric(scale(LogLength)),
    StdWeight = as.numeric(scale(LogWeight))) %>%
  select(Treatment3, Species_name, LogWidth, LogLength, LogWeight,
         StdWidth, StdLength, StdWeight)

seeds_plot <- seeds %>%
  group_by(plant_species, Treatment3) %>%
  filter(!plant_species %in% c("Theobroma_cacao", "Manihot_esculenta", "Artocarpus_heterophyllus", "Psidium_guajava",
                              "Borojoa_sp.", "Persea_americana")) %>%
  summarise(Treatment3 = first(Treatment3),
            Plot_ID = first(Plot_ID), RegTime = first(RegTime),
            Species_name = first(plant_species),
            .groups = "drop") 

### Seedlings
# BETTER CHECK SEED DISPERSAL SYNDROME OF SEEDLINGS

colSums(is.na(traits_seedlings)) # old has fewer NAs

seedlings_hbt <- seedlings %>%
  # filter(Plot_ID %in% master$Plot_ID[master$PREX == "PREX"]) %>%
  group_by(Species_name, Treatment3) %>%
  filter(!Species_name %in% c("Theobroma_cacao", "Manihot_esculenta", "Artocarpus_heterophyllus", "Psidium_guajava",
                               "Borojoa_sp.", "Persea_americana")) %>%
  filter(!Species_name %in% c("Coussapoa_villosa", "")) %>% # filter out non-zoochoric species based on Landim et al. 2026
  summarise(Treatment3 = first(Treatment3), 
            Species_name = first(Species_name), 
            .groups = 'drop') %>%
  left_join(traits_seedlings, by = "Species_name", relationship = "many-to-many") %>%
  ungroup()%>%
  mutate(Genus = word(Species_name, 1, sep = "_")) %>%
  group_by(Genus) %>%
  mutate(toughness = replace_na(data = toughness_old, replace = mean(toughness_old, na.rm = T)),
         thickness = replace_na(data = thickness_old, replace = mean(thickness_old, na.rm = T)),
         SLA = replace_na(data = SLA_old, replace = mean(SLA_old, na.rm = T)),
         LDMC = replace_na(data = LDMC_old, replace = mean(LDMC_old, na.rm = T))) %>%
  ungroup() %>%
  mutate(LogThough = log(toughness), LogThick = log(thickness)) %>%
  mutate(    
    StdThough = as.numeric(scale(LogThough)),
    StdThick = as.numeric(scale(LogThick)),
    StdSLA = as.numeric(scale(SLA)),
    StdLDMC = as.numeric(scale(LDMC))) %>%
  select(Treatment3, Species_name, LogThough, LogThick, SLA, LDMC,
         StdThough, StdThick, StdSLA, StdLDMC)

seedlings_plot <- seedlings %>%
  # filter(Plot_ID %in% master$Plot_ID[master$PREX == "PREX"]) %>%
  group_by(Species_name, Treatment3) %>%
  filter(!Species_name %in% c("Theobroma_cacao", "Manihot_esculenta", "Artocarpus_heterophyllus", "Psidium_guajava",
                               "Borojoa_sp.", "Persea_americana")) %>%
  summarise(Treatment3 = first(Treatment3),
            Plot_ID = first(Plot_ID), RegTime = first(RegTime),
            Species_name = first(Species_name),
            .groups = "drop") 

#### Functional spaces per habitat ####

#### Pollination

# Number of components to be used:
# cor <- cor.smooth(int_pol_hbt[, c("LogPL", "LogWL", "LogH", "LogCR")])
# eigen <- eigen(cor)
# permuted <- matrix(nrow=1000, ncol=4)
# for(i in 1:1000){
#   permuted_data <- apply(int_pol_hbt[, c("LogPL", "LogWL", "LogH", "LogCR")],2,sample)
#   permuted[i,] <- eigen(cor.smooth(permuted_data))$values
# }
# thresholds <- apply(permuted, 2, function(x) quantile(x, 0.95))
# print(eigen$values)
# print(thresholds)

pca_pol <- principal(int_pol_hbt[, c("StdPL", "StdWL", "StdH", "StdCR")], nfactor = 2, scores = TRUE, rotate = "varimax", covar = FALSE, missing = TRUE, use = "pairwise")

scores_pol <- as.data.frame(pca_pol$scores)
scores_pol$Treatment3 <- int_pol_hbt$Treatment3
scores_pol$interaction <- int_pol_hbt$interaction
scores_pol$group <- int_pol_hbt$group
scores_pol <- scores_pol %>% dplyr::select(group, interaction, Treatment3, RC1, RC2)

scores_pol$Treatment3 <-  factor(scores_pol$Treatment3,
                                 levels = c('old-growth forest', 'regeneration late', 'regeneration early'))

#### Seed-dispersal

# Number of components to be used:
cor <- cor.smooth(int_sd_hbt[, c("LogBM", "LogCM","LogGW", "LogFW", "HWIndex", "LogH")])
eigen <- eigen(cor)
permuted <- matrix(nrow=1000, ncol=6)
for(i in 1:1000){
  permuted_data <- apply(int_sd_hbt[, c("LogBM", "LogCM", "LogGW", "LogFW", "HWIndex", "LogH")],2,sample)
  permuted[i,] <- eigen(cor.smooth(permuted_data))$values
}
thresholds <- apply(permuted, 2, function(x) quantile(x, 0.95))
print(eigen$values)
print(thresholds)

pca_sd <- principal(int_sd_hbt[, c("StdBM", "StdCM", "StdGW", "StdFW", "StdHWI", "StdH")], nfactor = 2, scores = TRUE, rotate = "varimax", covar = FALSE, missing = TRUE, use = "pairwise")

scores_sd <- as.data.frame(pca_sd$scores)
scores_sd$Treatment3 <- int_sd_hbt$Treatment3
scores_sd$interaction <- int_sd_hbt$interaction
scores_sd$group <- int_sd_hbt$group
scores_sd <- scores_sd %>% dplyr::select(group, interaction, Treatment3, RC1, RC2)

scores_sd$Treatment3 <-  factor(scores_sd$Treatment3,
                                levels = c('old-growth forest', 'regeneration late', 'regeneration early'))

# pca_processes <- list(pca_pol, pca_sd)
# scores_processes <- list(scores_pol, scores_sd)
# saveRDS(pca_processes, file = "output/pca_processes.RDS")
# saveRDS(scores_processes, file = "output/scores_processes.RDS")

#### Seeds

# Number of components to be used:
# cor <- cor.smooth(seeds_hbt[, c("LogWidth", "LogLength", "LogWeight")])
# eigen <- eigen(cor)
# permuted <- matrix(nrow=1000, ncol=3)
# for(i in 1:1000){
#   permuted_data <- apply(seeds_hbt[, c("LogWidth", "LogLength", "LogWeight")],2,sample)
#   permuted[i,] <- eigen(cor.smooth(permuted_data))$values
# }
# thresholds <- apply(permuted, 2, function(x) quantile(x, 0.95))
# print(eigen$values)
# print(thresholds)

pca_seeds <- principal(seeds_hbt[, c("StdWidth", "StdLength", "StdWeight")], nfactor = 2, scores = TRUE, rotate = "varimax", covar = FALSE, missing = TRUE, use = "pairwise")

scores_seeds <- as.data.frame(pca_seeds$scores)
scores_seeds$Treatment3 <- seeds_hbt$Treatment3
scores_seeds$Species_name <- seeds_hbt$Species_name
scores_seeds <- scores_seeds %>% dplyr::select(Species_name, Treatment3, RC1, RC2)

scores_seeds$Treatment3 <-  factor(scores_seeds$Treatment3,
                                   levels = c('old-growth forest', 'regeneration late', 'regeneration early'))

#### Seedlings

# Number of components to be used:
# cor <- cor.smooth(seedlings_hbt[, c("LogThough", "LogThick", "SLA", "LDMC")])
# eigen <- eigen(cor)
# permuted <- matrix(nrow=1000, ncol=4)
# for(i in 1:1000){
#   permuted_data <- apply(seedlings_hbt[, c("LogThough", "LogThick", "SLA", "LDMC")],2,sample)
#   permuted[i,] <- eigen(cor.smooth(permuted_data))$values
# }
# thresholds <- apply(permuted, 2, function(x) quantile(x, 0.95))
# print(eigen$values)
# print(thresholds)

pca_sdlng <- principal(seedlings_hbt[, c("StdThough", "StdThick", "StdSLA", "StdLDMC")], nfactor = 2, scores = TRUE, rotate = "varimax", covar = FALSE, missing = TRUE, use = "pairwise")

scores_sdlng <- as.data.frame(pca_sdlng$scores)
scores_sdlng$Treatment3 <- seedlings_hbt$Treatment3
scores_sdlng$Species_name <- seedlings_hbt$Species_name
scores_sdlng <- scores_sdlng %>% dplyr::select(Species_name, Treatment3, RC1, RC2)

scores_sdlng$Treatment3 <-  factor(scores_sdlng$Treatment3,
                                levels = c('old-growth forest', 'regeneration late', 'regeneration early'))

#### Functional diversity ####

#### Pollination

orig_bees <- originality(int_bees_plot, scores_pol[scores_pol$group == "Bees", -1])
orig_moths <- originality(int_moths_plot, scores_pol[scores_pol$group == "Moths", -1])
orig_bat_pol <- originality(int_bat_pol_plot, scores_pol[scores_pol$group == "Bat_pol", -1])

FD_bees <- orig_bees %>%
  group_by(Plot_ID) %>%
  summarise(
    FDBees = mean(orig),
    FC1Bees = mean(RC1),
    FC2Bees = mean(RC2),
    .groups = "drop"
  ) 

FD_moths <- orig_moths %>%
  group_by(Plot_ID) %>%
  summarise(
    FDMoths = mean(orig),
    FC1Moths = mean(RC1),
    FC2Moths = mean(RC2),
    .groups = "drop"
  ) 

FD_bat_pol <- orig_bat_pol %>%
  group_by(Plot_ID) %>%
  summarise(
    FDBat_pol = mean(orig),
    FC1Bat_pol = mean(RC1),
    FC2Bat_pol = mean(RC2),
    .groups = "drop"
  ) 

#### Seed dispersal

orig_bats <- originality(int_bats_plot, scores_sd[scores_sd$group == "Bats", -1])
orig_birds <- originality(int_birds_plot, scores_sd[scores_sd$group == "Birds", -1])
orig_nf <- originality(int_nf_plot, scores_sd[scores_sd$group == "NF", -1])

FD_bats <- orig_bats %>%
  group_by(Plot_ID) %>%
  summarise(
    FDBats = mean(orig),
    FC1Bats = mean(RC1),
    FC2Bats = mean(RC2),
    .groups = "drop"
  ) 

FD_birds <- orig_birds %>%
  group_by(Plot_ID) %>%
  summarise(
    FDBirds = mean(orig),
    FC1Birds = mean(RC1),
    FC2Birds = mean(RC2),
    .groups = "drop"
  ) 

FD_nf <- orig_nf %>%
  group_by(Plot_ID) %>%
  summarise(
    FDNf = mean(orig),
    FC1Nf = mean(RC1),
    FC2Nf = mean(RC2),
    .groups = "drop"
  ) 

#### Seeds

orig_seeds <- originality(seeds_plot, scores_seeds)

TD_seeds <- seeds %>%
  rename(Species_name = plant_species) %>%
  group_by(Plot_ID) %>%
  summarise(AbSeeds = n(),
            RichSeeds = n_distinct(Species_name),
            .groups = "drop")
Shn_seeds <- seeds %>%
  rename(Species_name = plant_species) %>%
  count(Plot_ID, Species_name) %>%
  group_by(Plot_ID) %>%
  summarise(ShnSeeds = diversity(n, index = "shannon"),
            .groups = "drop")


FD_seeds <- orig_seeds %>%
  group_by(Plot_ID) %>%
  summarise(
    FDSeeds = mean(orig),
    FC1Seeds = mean(RC1),
    FC2Seeds = mean(RC2),
    .groups = "drop"
  ) 


#### Seedlings

orig_sdlngs <- originality(seedlings_plot, scores_sdlng)

TD_sdlngs <- seedlings %>%
  group_by(Plot_ID) %>%
  summarise(AbSdlng = n(),
            RichSdlng = n_distinct(Species_name),
            .groups = "drop")
Shn_sdlngs <- seedlings %>%
  count(Plot_ID, Species_name) %>%
  group_by(Plot_ID) %>%
  summarise(ShnSdlng = diversity(n, index = "shannon"),
            .groups = "drop")


FD_sdlngs <- orig_sdlngs %>%
  group_by(Plot_ID) %>%
  summarise(
    FDSdlng = mean(orig),
    FC1Sdlng = mean(RC1),
    FC2Sdlng = mean(RC2),
    .groups = "drop"
  ) 


#### Explanatory variables (recovery time, forest structure and connectivity) ####

expl <- read.csv(file = "data/processed/env_variables.csv")

## Connectivity

expl_con <- expl %>% select(Treatment3, Plot_ID, RegTime, Forest_1km, Forest_500m, Forest_100m, Distance_forest, Distance_edge, Cacao_1km, Pasture_1km, Cacao_Reg1_1km, Pasture_Reg1_1km, Cacao_Reg2_1km, Pasture_Reg2_1km)
expl_con$Distance_edge <- -expl_con$Distance_edge
expl_con$Distance <- expl_con$Distance_forest + expl_con$Distance_edge

# hist(exp(expl_con$Forest_1km))
# hist(exp(expl_con$Forest_500m))
# hist(log1p(expl_con$Forest_100m))
# hist((expl_con$Distance+max(abs(expl_con$Distance))))
# hist(log1p(expl_con$Cacao_1km))
# hist(log1p(expl_con$Pasture_1km))
# hist(log1p(expl_con$Cacao_Reg1_1km))
# hist(log1p(expl_con$Pasture_Reg1_1km))
# hist(log1p(expl_con$Cacao_Reg2_1km))
# hist(log1p(expl_con$Pasture_Reg2_1km))

con <- expl_con %>%
  mutate(
    t_100m = (Forest_100m),
    t_500m = exp(Forest_500m),
    t_1km = exp(Forest_1km),
    t_dist = (Distance+max(abs(Distance))),
    t_cacao = log1p(Cacao_1km),
    t_pasture = log1p(Pasture_1km),
    t_cacao1 = log1p(Cacao_Reg1_1km),
    t_pasture1 = log1p(Pasture_Reg1_1km),
    t_cacao2 = log1p(Cacao_Reg2_1km),
    t_pasture2 = log1p(Pasture_Reg2_1km),
    Std_100m = scale((Forest_100m)),
    Std_500m = scale(exp(Forest_500m)),
    Std_1km = scale(exp(Forest_1km)),
    Std_dist = scale((Distance+max(abs(Distance)))),
    Std_cacao = scale(log1p(Cacao_1km)),
    Std_pasture = scale(log1p(Pasture_1km)),
    Std_cacao1 = scale(log1p(Cacao_Reg1_1km)),
    Std_pasture1 = scale(log1p(Pasture_Reg1_1km)),
    Std_cacao2 = scale(log1p(Cacao_Reg2_1km)),
    Std_pasture2 = scale(log1p(Pasture_Reg2_1km)),
  )

# Number of components to be used:
# cor <- cor.smooth(con[, c("t_1km", "t_500m", "t_100m", "t_dist", "t_cacao", "t_pasture", "t_cacao1", "t_pasture1", "t_cacao2", "t_pasture2")])
# eigen <- eigen(cor)
# permuted <- matrix(nrow=1000, ncol=10)
# for(i in 1:1000){
#   permuted_data <- apply(con[, c("t_1km", "t_500m", "t_100m", "t_dist", "t_cacao", "t_pasture", "t_cacao1", "t_pasture1", "t_cacao2", "t_pasture2")],2,sample)
#   permuted[i,] <- eigen(cor.smooth(permuted_data))$values
# }
# thresholds <- apply(permuted, 2, function(x) quantile(x, 0.95))
# print(eigen$values)
# print(thresholds)

con_pca <- principal(con[, c("Std_100m", "Std_500m", "Std_1km", "Std_dist", "Std_cacao", "Std_pasture", "Std_cacao1", "Std_pasture1", "Std_cacao2", "Std_pasture2")], nfactor = 1, scores = TRUE, rotate = "varimax", covar = FALSE, missing = TRUE, use = "pairwise")

con_idx <- as.data.frame(con_pca$scores*-1) 
con_idx$Plot_ID <- expl_con$Plot_ID

con_idx <- con_idx %>% rename(ConIndex = PC1) %>% select(Plot_ID, ConIndex)

## Forest structure

expl_str <- expl %>% select(Plot_ID, PREX, Treatment3, VerticalVH, Max_tree_height, AGB_wild)
str(expl_str)

# hist(exp(expl_str$VerticalVH))
# hist((expl_str$Max_tree_height))
# hist(log1p(expl_str$AGB_wild))

strc <- expl_str %>%
  mutate(
    t_vvh = exp(VerticalVH),
    t_agb = log1p(AGB_wild),
    Std_vvh = scale(exp(VerticalVH)),
    Std_mth = scale((Max_tree_height)),
    Std_agb = scale(log1p(AGB_wild))
  )

# Number of components to be used:
# cor <- cor.smooth(strc[, c("t_vvh", "Max_tree_height", "t_agb")])
# eigen <- eigen(cor)
# permuted <- matrix(nrow=1000, ncol=3)
# for(i in 1:1000){
#   permuted_data <- apply(strc[, c("t_vvh", "Max_tree_height", "t_agb")],2,sample)
#   permuted[i,] <- eigen(cor.smooth(permuted_data))$values
# }
# thresholds <- apply(permuted, 2, function(x) quantile(x, 0.95))
# print(eigen$values)
# print(thresholds)

str_pca <- principal(strc[, c("Std_vvh", "Std_mth", "Std_agb")], nfactor = 1, scores = TRUE, rotate = "varimax", covar = FALSE, missing = TRUE, use = "pairwise")

str_idx <- as.data.frame(str_pca$scores)
str_idx$Plot_ID <- strc$Plot_ID

str_idx <- str_idx %>% rename(StrIndex = PC1) %>% select(Plot_ID, StrIndex)


#### Final data frame ####

model_df <- expl %>% select(Treatment3, Plot_ID, RegTime) %>%
  left_join(con_idx, by = "Plot_ID") %>%
  left_join(str_idx, by = "Plot_ID") %>%
  left_join(FD_bees, by = "Plot_ID") %>%
  left_join(FD_moths, by = "Plot_ID") %>%
  left_join(FD_bat_pol, by = "Plot_ID") %>%
  left_join(FD_bats, by = "Plot_ID") %>%
  left_join(FD_birds, by = "Plot_ID") %>%
  left_join(FD_nf, by = "Plot_ID")# %>%
  # left_join(FD_seeds, by = "Plot_ID") %>%
  # left_join(TD_seeds, by = "Plot_ID") %>%
  # left_join(Shn_seeds, by = "Plot_ID") %>%
  # left_join(FD_sdlngs, by = "Plot_ID") %>%
  # left_join(TD_sdlngs, by = "Plot_ID") %>%
  # left_join(Shn_sdlngs, by = "Plot_ID") %>%
  # mutate(
  #   across(
  #     c(AbSdlng, RichSdlng, ShnSdlng),
  #     ~ tidyr::replace_na(.x, 0)
  #   )
  # )

str(model_df)


# saveRDS(list(scores_pol, scores_sd), file = here::here("output", "scores_processes.3.RDS"))
# saveRDS(list(pca_pol, pca_sd), file = here::here("output", "pca_processes.3.RDS"))

# write.csv(model_df, file = here::here("data", "processed", "model_df.3.csv"), row.names = F)
