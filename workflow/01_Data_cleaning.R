#### Cleaning data sets ####

library(dplyr)
library(tidyr)

#### Interactions ####

master <- read.csv(file = "data/raw/Ecuador Plots_Master_V2.csv")
master$RegTime <- 2023 - master$Regeneration_year
master$RegTime[master$Treatment == "active cacao"] <- 0
master$RegTime[master$Treatment == "active pasture"] <- 0
master$RegTime[master$Treatment == "old-growth forest"] <- 55

regtime <- master %>% select(Plot_ID, RegTime)
  
int_do <- read.csv(file = "data/raw/SP4/int_direct.obs_org.csv")
int_do$Plot_ID <- as.factor(int_do$Plot_ID)
int_do$Treatment2 <- as.factor(int_do$Treatment2)
int_do <- int_do %>% select(Plot_ID, taxon, animal_species, plant_species, type, unit, RegTime, Treatment2, method)

int_ct <- read.csv(file = "data/raw/SP4/int_camera.trap_org.csv")
int_ct$Plot_ID <- as.factor(int_ct$Plot_ID)
int_ct$Treatment2 <- as.factor(int_ct$Treatment2)
int_ct <- int_ct %>% select(Plot_ID, taxon, animal_species, plant_species, type, unit, RegTime, Treatment2, method)

int_bat <- read.csv(file = "data/raw/SP4/Bat - Seed_ Interaction_SE.csv")
int_bat$Plot_ID <- as.factor(int_bat$Plot_ID)
int_bat$Legacy_7G <- as.factor(int_bat$Legacy_7G)
int_bat <- int_bat %>%  
  left_join(regtime, by = "Plot_ID") %>%
  rename(animal_species = Bat_sp, plant_species = Seed_sp, Treatment2 = Legacy_7G) %>%
  mutate(taxon  = "Chiroptera", type   = "feces", unit   = "seed", method = "mist_net") %>%
  select(Plot_ID, taxon, animal_species, plant_species, Seed_Morpho, type, unit, RegTime, Treatment2, method)

int_bat[int_bat$animal_species == "Sturnia_luisi", "animal_species"] <- "Sturnira_luisi"

freq_bat <- read.csv(file = "data/raw/SP4/Bats_plots_ SE.csv")
freq_bat_long <- freq_bat %>%
  pivot_longer(
    -Plot_ID,
    names_to  = "animal_species",
    values_to = "n"
  ) %>%
  filter(!is.na(n), n > 0)

int_bat <- int_bat %>%
  left_join(freq_bat_long, by = c("Plot_ID", "animal_species")) %>%
  uncount(weights = n)

#### Plant names ####

Plants <- read.csv(file = "data/raw/SP4/traits_plants_org.csv")[-c(1:3),-1]
Plants[Plants$species=="Nectandra_purpurea_cf.", "species"] <- "Nectandra_purpurea"
Plants[Plants$species=="Osteophloeum_platis_cf.", "species"] <- "Osteophloeum_platyspermum"
Plants[Plants$species=="Miconia_'alargada tres venas abajo'", "species"] <- "Miconia_multiplicata"
Plants[Plants$species=="Miconia 'alargada tres venas abajo'", "species"] <- "Miconia_multiplicata"
Plants[Plants$species=="A-PA5603", "species"] <- "Lantana_sp."
Plants[Plants$species=="Guarea guidonia", "species"] <- "Guarea_guidonia"
Plants[Plants$species=="Myrcia_cf._aliena", "species"] <- "Myrcia_aliena"
Plants[Plants$species=="Ficus_cf._pertusa", "species"] <- "Ficus_pertusa"
Plants[Plants$species=="Ficus_'matapalo'", "species"] <- "Ficus_brevibracteata"
Plants[Plants$species=="dormilon", "species"] <- "Hymenolobium_heterocarpum"
Plants[Plants$species=="Miconia_'peciolo rojo'", "species"] <- "Miconia_multiplicata"
Plants[Plants$species=="Plinia_'longifolia'", "species"] <- "Eugenia_multirimosa"
Plants[Plants$species=="OG45-71", "species"] <- "Couepia_subcordata"
Plants[Plants$species=="OG37-176", "species"] <- "Eschweilera_decolorans"
Plants[Plants$species=="Conostegia_'rosada'", "species"] <- "Miconia_conocuatrecasii"
Plants[Plants$species=="Meliosma_cf._herbertii", "species"] <- "Meliosma_herbertii"
Plants[Plants$species=="Miconia_cf._montana", "species"] <- "Miconia_montana"
Plants[Plants$species=="Ficus_sp.", "species"] <- "Ficus_maxima"
Plants[Plants$species=="OG41-25", "species"] <- "Caryocar_glabrum"
Plants[Plants$species=="OG38-85", "species"] <- "Coussarea_latifolia"
Plants[Plants$species=="Borojoa sp.", "species"] <- "Borojoa_sp."
Plants[Plants$species=="Nectandra_purpurea", "species"] <- "Damburneya_purpurea"
Plants[Plants$species=="Cordia_sp.", "species"] <- "Exarata_chocoensis"
Plants[Plants$species=="Miconia_oraria", "species"] <- "Miconia_conocuatrecasii"
Plants[Plants$species=="Conostegia_superba", "species"] <- "Miconia_'grande'"
Plants[Plants$species=="Licania_cf._dodsonii", "species"] <- "Licania_glauca"
Plants[Plants$species=="Guarea_cf._bullata_cf._kunthiana", "species"] <- "Guarea_macrophylla"
Plants[Plants$species=="OG41-77", "species"] <- "Turpinia_occidentalis"
Plants[Plants$species=="Guarea_cf._cartaguenya", "species"] <- "Guarea_cartaguenya"
Plants[Plants$species=="Croton tessmannii", "species"] <- "Croton_tessmannii"
Plants[Plants$species=="crasa'", "species"] <- "Daphnopsis_occulta"
Plants[Plants$species=="Gustavia_aff._hexapetala", "species"] <- "Gustavia_hexapetala"
Plants[Plants$species=="Eschweilera???", "species"] <- "Klarobelia_megalocarpa"
Plants[Plants$species=="Piperacea_sp4", "species"] <- "Piper_grande"
Plants[Plants$species=="A-OG49002", "species"] <- "Henriettella_lawrancei"
Plants[Plants$species=="Piperacea_sp5", "species"] <- "Piper_pseudonobile"
Plants[Plants$species=="Piperacea_sp2", "species"] <- "Piper_pseudonobile"
Plants[Plants$species=="Miconia_'peluda'", "species"] <- "Miconia_cooperi"
Plants[Plants$species=="Miconia_pequena", "species"] <- "Miconia_cooperi"
Plants[Plants$species=="Tabebuia_rosea", "species"] <- "Handroanthus_chrysanthus"
Plants[Plants$species=="Cordia_cf._'achiote'", "species"] <- "Bixa_orellana"
Plants[Plants$species=="Piper_pseudonobile", "species"] <- "Piper_pequeina"
Plants[Plants$species=="Piperacea_sp3", "species"] <- "Piper_pequeina"


int_do[int_do$plant_species=="Nectandra_purpurea_cf.", "plant_species"] <- "Nectandra_purpurea"
int_do[int_do$plant_species=="Osteophloeum_platis_cf.", "plant_species"] <- "Osteophloeum_platyspermum"
int_do[int_do$plant_species=="Miconia_'alargada tres venas abajo'", "plant_species"] <- "Miconia_multiplicata"
int_do[int_do$plant_species=="Miconia 'alargada tres venas abajo'", "plant_species"] <- "Miconia_multiplicata"
int_do[int_do$plant_species=="A-PA5603", "plant_species"] <- "Lantana_sp."
int_do[int_do$plant_species=="Guarea guidonia", "plant_species"] <- "Guarea_guidonia"
int_do[int_do$plant_species=="Myrcia_cf._aliena", "plant_species"] <- "Myrcia_aliena"
int_do[int_do$plant_species=="Ficus_cf._pertusa", "plant_species"] <- "Ficus_pertusa"
int_do[int_do$plant_species=="Ficus_'matapalo'", "plant_species"] <- "Ficus_brevibracteata"
int_do[int_do$plant_species=="dormilon", "plant_species"] <- "Hymenolobium_heterocarpum"
int_do[int_do$plant_species=="Miconia_'peciolo rojo'", "plant_species"] <- "Miconia_multiplicata"
int_do[int_do$plant_species=="Plinia_'longifolia'", "plant_species"] <- "Eugenia_multirimosa"
int_do[int_do$plant_species=="OG45-71", "plant_species"] <- "Couepia_subcordata"
int_do[int_do$plant_species=="OG37-176", "plant_species"] <- "Eschweilera_decolorans"
int_do[int_do$plant_species=="Conostegia_'rosada'", "plant_species"] <- "Miconia_conocuatrecasii"
int_do[int_do$plant_species=="Meliosma_cf._herbertii", "plant_species"] <- "Meliosma_herbertii"
int_do[int_do$plant_species=="Miconia_cf._montana", "plant_species"] <- "Miconia_montana"
int_do[int_do$plant_species=="Ficus_sp.", "plant_species"] <- "Ficus_maxima"
int_do[int_do$plant_species=="OG41-25", "plant_species"] <- "Caryocar_glabrum"
int_do[int_do$plant_species=="OG38-85", "plant_species"] <- "Coussarea_latifolia"
int_do[int_do$plant_species=="Borojoa sp.", "plant_species"] <- "Borojoa_sp."
int_do[int_do$plant_species=="Nectandra_purpurea", "plant_species"] <- "Damburneya_purpurea"
int_do[int_do$plant_species=="Cordia_sp.", "plant_species"] <- "Exarata_chocoensis"
int_do[int_do$plant_species=="Miconia_oraria", "plant_species"] <- "Miconia_conocuatrecasii"
int_do[int_do$plant_species=="Conostegia_superba", "plant_species"] <- "Miconia_'grande'"
int_do[int_do$plant_species=="Licania_cf._dodsonii", "plant_species"] <- "Licania_glauca"
int_do[int_do$plant_species=="Guarea_cf._bullata_cf._kunthiana", "plant_species"] <- "Guarea_macrophylla"
int_do[int_do$plant_species=="OG41-77", "plant_species"] <- "Turpinia_occidentalis"
int_do[int_do$plant_species=="Guarea_cf._cartaguenya", "plant_species"] <- "Guarea_cartaguenya"
int_do[int_do$plant_species=="Croton tessmannii", "plant_species"] <- "Croton_tessmannii"
int_do[int_do$plant_species=="crasa'", "plant_species"] <- "Daphnopsis_occulta"
int_do[int_do$plant_species=="Gustavia_aff._hexapetala", "plant_species"] <- "Gustavia_hexapetala"
int_do[int_do$plant_species=="Eschweilera???", "plant_species"] <- "Klarobelia_megalocarpa"
int_do[int_do$plant_species=="Piperacea_sp4", "plant_species"] <- "Piper_grande"
int_do[int_do$plant_species=="A-OG49002", "plant_species"] <- "Henriettella_lawrancei"
int_do[int_do$plant_species=="Piperacea_sp5", "plant_species"] <- "Piper_pseudonobile"
int_do[int_do$plant_species=="Piperacea_sp2", "plant_species"] <- "Piper_pseudonobile"
int_do[int_do$plant_species=="Miconia_'peluda'", "plant_species"] <- "Miconia_cooperi"
int_do[int_do$plant_species=="Miconia_pequena", "plant_species"] <- "Miconia_cooperi"
int_do[int_do$plant_species=="Tabebuia_rosea", "plant_species"] <- "Handroanthus_chrysanthus"

int_ct[int_ct$plant_species=="Nectandra_purpurea_cf.", "plant_species"] <- "Nectandra_purpurea"
int_ct[int_ct$plant_species=="Osteophloeum_platis_cf.", "plant_species"] <- "Osteophloeum_platyspermum"
int_ct[int_ct$plant_species=="Miconia_'alargada tres venas abajo'", "plant_species"] <- "Miconia_multiplicata"
int_ct[int_ct$plant_species=="Miconia 'alargada tres venas abajo'", "plant_species"] <- "Miconia_multiplicata"
int_ct[int_ct$plant_species=="A-PA5603", "plant_species"] <- "Lantana_sp."
int_ct[int_ct$plant_species=="Guarea guidonia", "plant_species"] <- "Guarea_guidonia"
int_ct[int_ct$plant_species=="Myrcia_cf._aliena", "plant_species"] <- "Myrcia_aliena"
int_ct[int_ct$plant_species=="Ficus_cf._pertusa", "plant_species"] <- "Ficus_pertusa"
int_ct[int_ct$plant_species=="Ficus_'matapalo'", "plant_species"] <- "Ficus_brevibracteata"
int_ct[int_ct$plant_species=="dormilon", "plant_species"] <- "Hymenolobium_heterocarpum"
int_ct[int_ct$plant_species=="Miconia_'peciolo rojo'", "plant_species"] <- "Miconia_multiplicata"
int_ct[int_ct$plant_species=="Plinia_'longifolia'", "plant_species"] <- "Eugenia_multirimosa"
int_ct[int_ct$plant_species=="OG45-71", "plant_species"] <- "Couepia_subcordata"
int_ct[int_ct$plant_species=="OG37-176", "plant_species"] <- "Eschweilera_decolorans"
int_ct[int_ct$plant_species=="Conostegia_'rosada'", "plant_species"] <- "Miconia_conocuatrecasii"
int_ct[int_ct$plant_species=="Meliosma_cf._herbertii", "plant_species"] <- "Meliosma_herbertii"
int_ct[int_ct$plant_species=="Miconia_cf._montana", "plant_species"] <- "Miconia_montana"
int_ct[int_ct$plant_species=="Ficus_sp.", "plant_species"] <- "Ficus_maxima"
int_ct[int_ct$plant_species=="OG41-25", "plant_species"] <- "Caryocar_glabrum"
int_ct[int_ct$plant_species=="OG38-85", "plant_species"] <- "Coussarea_latifolia"
int_ct[int_ct$plant_species=="Borojoa sp.", "plant_species"] <- "Borojoa_sp."
int_ct[int_ct$plant_species=="Nectandra_purpurea", "plant_species"] <- "Damburneya_purpurea"
int_ct[int_ct$plant_species=="Cordia_sp.", "plant_species"] <- "Exarata_chocoensis"
int_ct[int_ct$plant_species=="Miconia_oraria", "plant_species"] <- "Miconia_conocuatrecasii"
int_ct[int_ct$plant_species=="Conostegia_superba", "plant_species"] <- "Miconia_'grande'"
int_ct[int_ct$plant_species=="Licania_cf._dodsonii", "plant_species"] <- "Licania_glauca"
int_ct[int_ct$plant_species=="Guarea_cf._bullata_cf._kunthiana", "plant_species"] <- "Guarea_macrophylla"
int_ct[int_ct$plant_species=="OG41-77", "plant_species"] <- "Turpinia_occidentalis"
int_ct[int_ct$plant_species=="Guarea_cf._cartaguenya", "plant_species"] <- "Guarea_cartaguenya"
int_ct[int_ct$plant_species=="Croton tessmannii", "plant_species"] <- "Croton_tessmannii"
int_ct[int_ct$plant_species=="crasa'", "plant_species"] <- "Daphnopsis_occulta"
int_ct[int_ct$plant_species=="Gustavia_aff._hexapetala", "plant_species"] <- "Gustavia_hexapetala"
int_ct[int_ct$plant_species=="Eschweilera???", "plant_species"] <- "Klarobelia_megalocarpa"
int_ct[int_ct$plant_species=="Piperacea_sp4", "plant_species"] <- "Piper_grande"
int_ct[int_ct$plant_species=="A-OG49002", "plant_species"] <- "Henriettella_lawrancei"
int_ct[int_ct$plant_species=="Piperacea_sp5", "plant_species"] <- "Piper_pseudonobile"
int_ct[int_ct$plant_species=="Piperacea_sp2", "plant_species"] <- "Piper_pseudonobile"
int_ct[int_ct$plant_species=="Miconia_'peluda'", "plant_species"] <- "Miconia_cooperi"
int_ct[int_ct$plant_species=="Miconia_pequena", "plant_species"] <- "Miconia_cooperi"
int_ct[int_ct$plant_species=="Tabebuia_rosea", "plant_species"] <- "Handroanthus_chrysanthus"

int_bat[int_bat$Seed_Morpho=="M_06", "plant_species"] <- "Piptocoma_discolor"
int_bat[int_bat$Seed_Morpho=="M_08", "plant_species"] <- "Fusispermum_minutiflorum"
int_bat[int_bat$Seed_Morpho=="M_10", "plant_species"] <- "Neosprucea_pedicellata"
int_bat[int_bat$Seed_Morpho=="M_100", "plant_species"] <- "Byrsonima_ligustrifolia"
int_bat[int_bat$Seed_Morpho=="M_101", "plant_species"] <- "Piper_grande"
int_bat[int_bat$Seed_Morpho=="M_102", "plant_species"] <- "Faramea_langlassei"
int_bat[int_bat$Seed_Morpho=="M_104", "plant_species"] <- "Nectandra_paucinervia"
int_bat[int_bat$Seed_Morpho=="M_105", "plant_species"] <- "Piper_pequeina"
int_bat[int_bat$Seed_Morpho=="M_106", "plant_species"] <- "Beilschmiedia_costaricensis"
int_bat[int_bat$Seed_Morpho=="M_108", "plant_species"] <- "Inga_oerstediana"
int_bat[int_bat$Seed_Morpho=="M_11", "plant_species"] <- "Ficus_pertusa"
int_bat[int_bat$Seed_Morpho=="M_110", "plant_species"] <- "Lozania_klugii"
int_bat[int_bat$Seed_Morpho=="M_111", "plant_species"] <- "Pouteria_leptopedicellata"
int_bat[int_bat$Seed_Morpho=="M_113", "plant_species"] <- "Byrsonima_ligustrifolia"
int_bat[int_bat$Seed_Morpho=="M_114", "plant_species"] <- "Piper_pequeina"
int_bat[int_bat$Seed_Morpho=="M_12", "plant_species"] <- "Ficus_brevibracteata"
int_bat[int_bat$Seed_Morpho=="M_13", "plant_species"] <- "Piptocoma_discolor"
int_bat[int_bat$Seed_Morpho=="M_14", "plant_species"] <- "Piptocoma_discolor"
int_bat[int_bat$Seed_Morpho=="M_15", "plant_species"] <- "Banara_guianensis"
int_bat[int_bat$Seed_Morpho=="M_16", "plant_species"] <- "Banara_guianensis"
int_bat[int_bat$Seed_Morpho=="M_19", "plant_species"] <- "Minquartia_guianensis"
int_bat[int_bat$Seed_Morpho=="M_20", "plant_species"] <- "Piper_pequeina"
int_bat[int_bat$Seed_Morpho=="M_21", "plant_species"] <- "Dendropanax_oblonga"
int_bat[int_bat$Seed_Morpho=="M_23", "plant_species"] <- "Dacryodes_cupularis"
int_bat[int_bat$Seed_Morpho=="M_24", "plant_species"] <- "Cordia_alliodora"
int_bat[int_bat$Seed_Morpho=="M_26", "plant_species"] <- "Exarata_chocoensis"
int_bat[int_bat$Seed_Morpho=="M_27", "plant_species"] <- "Naucleopsis_chiguila"
int_bat[int_bat$Seed_Morpho=="M_33", "plant_species"] <- "Cecropia_insignis"
int_bat[int_bat$Seed_Morpho=="M_36", "plant_species"] <- "Citronella_incarum"
int_bat[int_bat$Seed_Morpho=="M_37", "plant_species"] <- "Ficus_pertusa"
int_bat[int_bat$Seed_Morpho=="M_38", "plant_species"] <- "Discophora_guianensis"
int_bat[int_bat$Seed_Morpho=="M_39", "plant_species"] <- "Piper_pequeina"
int_bat[int_bat$Seed_Morpho=="M_42", "plant_species"] <- "Pouteria_gigantea"
int_bat[int_bat$Seed_Morpho=="M_44", "plant_species"] <- "Piper_grande"
int_bat[int_bat$Seed_Morpho=="M_45", "plant_species"] <- "Casearia_arborea"
int_bat[int_bat$Seed_Morpho=="M_46", "plant_species"] <- "Piptocoma_discolor"
int_bat[int_bat$Seed_Morpho=="M_49", "plant_species"] <- "Minquartia_guianensis"
int_bat[int_bat$Seed_Morpho=="M_50", "plant_species"] <- "Dendropanax_oblonga"
int_bat[int_bat$Seed_Morpho=="M_51", "plant_species"] <- "Piper_grande"
int_bat[int_bat$Seed_Morpho=="M_52", "plant_species"] <- "Cecropia_insignis"
int_bat[int_bat$Seed_Morpho=="M_54", "plant_species"] <- "Saurauia_herthae"
int_bat[int_bat$Seed_Morpho=="M_55", "plant_species"] <- "Piper_grande"
int_bat[int_bat$Seed_Morpho=="M_56", "plant_species"] <- "Aegiphila_alba"
int_bat[int_bat$Seed_Morpho=="M_60", "plant_species"] <- "Coussapoa_contorta"
int_bat[int_bat$Seed_Morpho=="M_61", "plant_species"] <- "Dendropanax_oblonga"
int_bat[int_bat$Seed_Morpho=="M_63", "plant_species"] <- "Miconia_dorsiloba"
int_bat[int_bat$Seed_Morpho=="M_64", "plant_species"] <- "Piper_pequeina"
int_bat[int_bat$Seed_Morpho=="M_66", "plant_species"] <- "Genipa_americana"
int_bat[int_bat$Seed_Morpho=="M_67", "plant_species"] <- "Nectandra_sp.2"
int_bat[int_bat$Seed_Morpho=="M_69", "plant_species"] <- "Bunchosia_cornifolia"
int_bat[int_bat$Seed_Morpho=="M_71", "plant_species"] <- "Trema_micrantha"
int_bat[int_bat$Seed_Morpho=="M_72", "plant_species"] <- "Banara_guianensis"
int_bat[int_bat$Seed_Morpho=="M_73", "plant_species"] <- "Apeiba_membranacea"
int_bat[int_bat$Seed_Morpho=="M_74", "plant_species"] <- "Psidium_guajava"
int_bat[int_bat$Seed_Morpho=="M_75", "plant_species"] <- "Psidium_guajava"
int_bat[int_bat$Seed_Morpho=="M_78", "plant_species"] <- "Hieronyma_macrocarpa"
int_bat[int_bat$Seed_Morpho=="M_79", "plant_species"] <- "Trema_integerrima"
int_bat[int_bat$Seed_Morpho=="M_83", "plant_species"] <- "Piper_pequeina"
int_bat[int_bat$Seed_Morpho=="M_85", "plant_species"] <- "Piper_grande"
int_bat[int_bat$Seed_Morpho=="M_87", "plant_species"] <- "Coussapoa_contorta"
int_bat[int_bat$Seed_Morpho=="M_89", "plant_species"] <- "Pouteria_caimito"
int_bat[int_bat$Seed_Morpho=="M_91", "plant_species"] <- "Piper_grande"
int_bat[int_bat$Seed_Morpho=="M_95", "plant_species"] <- "Trema_micrantha"
int_bat[int_bat$Seed_Morpho=="M_98", "plant_species"] <- "Miconia_multiplicata"
int_bat[int_bat$Seed_Morpho=="M_99", "plant_species"] <- "Piptocoma_discolor"
int_bat[int_bat$plant_species =="Bunchosia_cornifolia", "plant_species"] <- "Bunchosia_nitida"
int_bat[int_bat$plant_species =="Conostegia_cuatrecasasii", "plant_species"] <- "Miconia_conocuatrecasii"

#### Habitat types ####

int_do <- int_do %>%
  mutate(Treatment3 = as.factor(case_when(
    Treatment2 %in% c("active cacao", "active pasture") ~ "regeneration early",
    Treatment2 %in% c("cacao regeneration early", "pasture regeneration early") ~ "regeneration early",
    Treatment2 %in% c("cacao regeneration late", "pasture regeneration late") ~ "regeneration late",
    TRUE ~ as.character(Treatment2) 
  )))

int_ct <- int_ct %>%
  mutate(Treatment3 = as.factor(case_when(
    Treatment2 %in% c("active cacao", "active pasture") ~ "regeneration early",
    Treatment2 %in% c("cacao regeneration early", "pasture regeneration early") ~ "regeneration early",
    Treatment2 %in% c("cacao regeneration late", "pasture regeneration late") ~ "regeneration late",
    TRUE ~ as.character(Treatment2) 
  )))

int_bat <- int_bat %>%
  mutate(Treatment3 = as.factor(case_when(
    Treatment2 %in% c("Cacao", "Pasture") ~ "regeneration early",
    Treatment2 %in% c("CR_early", "PR_early") ~ "regeneration early",
    Treatment2 %in% c("CR_late", "PR_late") ~ "regeneration late",
    Treatment2 %in% "Old_growth" ~ "old-growth forest",
    TRUE ~ as.character(Treatment2) 
  )))

Plants <- Plants %>%
  mutate(Treatment3 = as.factor(case_when(
    Treatment2 %in% c("active cacao", "active pasture") ~ "regeneration early",
    Treatment2 %in% c("cacao regeneration early", "pasture regeneration early") ~ "regeneration early",
    Treatment2 %in% c("cacao regeneration late", "pasture regeneration late") ~ "regeneration late",
    TRUE ~ as.character(Treatment2) 
  )))

#### Organizing per functional group ####

do_ct <- bind_rows(int_do, int_ct)
int_birds <- do_ct %>% filter(taxon == "Birds") %>% distinct(Plot_ID, animal_species, plant_species, .keep_all = TRUE)
int_nf <- do_ct %>% filter(!taxon %in% c("Birds", "Chiroptera")) %>% distinct(Plot_ID, animal_species, plant_species, .keep_all = TRUE)
int_bat <- int_bat %>% distinct(Plot_ID, animal_species, plant_species, .keep_all = TRUE)

#### Traits ####

traits_do <- read.csv(file = "data/raw/SP4/traits_direct.obs_org.csv")
traits_do$Plot_ID <- as.factor(traits_do$Plot_ID)

traits_ct <- read.csv(file = "data/raw/SP4/traits_cam.trap_org.csv")
traits_ct$Plot_ID <- as.factor(traits_ct$Plot_ID)

traits_bat <- read.csv(file = "data/raw/SP4/Functional_Bat_matrix_SE.csv")
traits_bat <- traits_bat %>% rename(animal_species = specie, GapeWidth = JW, HWIndex = HWI, BodyMass = W) %>%
  select(animal_species, GapeWidth, HWIndex, BodyMass)

do_ct_traits <- bind_rows(
  traits_ct %>% mutate(source = "ct"),
  traits_do %>% mutate(source = "do"))

bird_trait_map <- tibble::tribble(
  ~trait,             ~new_col,
  "Beak.Width",       "BeakWidth",
  "Hand.wing.Index",  "HWIndex",
  "BodyMass",         "BodyMass")
traits_birds <- do_ct_traits %>%
  filter(taxon == "Birds") %>%
  inner_join(bird_trait_map, by = "trait") %>%
  group_by(species, new_col) %>%
  summarise(value = first(mean_value_numeric), .groups = "drop") %>%
  pivot_wider(
    names_from = new_col,
    values_from = value) %>%
  rename(animal_species = species)

nf_trait_map <- tibble::tribble(
  ~trait,          ~new_col,
  "mandibleWidth", "GapeWidth",
  "BodyMass",      "BodyMass")
traits_nf <- do_ct_traits %>%
  filter(!taxon %in% c("Birds", "Chiroptera")) %>%
  inner_join(nf_trait_map, by = "trait") %>%
  group_by(species, new_col) %>%
  summarise(value = first(mean_value_numeric), .groups = "drop") %>%
  pivot_wider(
    names_from = new_col,
    values_from = value) %>%
  rename(animal_species = species)

Plants <- Plants %>%
  select(species, trait, value_numeric)

bat_plants <- read.csv("data/raw/Bat_plants.csv")

bat_plants <- bat_plants %>% 
  transmute(plant_species, FruitWidth, Height, CropMass = NA_real_)

fruit_width <- Plants %>%
  filter(trait == "FruitWidth") %>%
  group_by(species) %>%
  summarize(FruitWidth = mean(value_numeric, na.rm = TRUE)) %>%
  dplyr::select(plant_species = species, FruitWidth) %>%
  distinct()

height <- Plants %>%
  filter(trait == "height") %>%
  group_by(species) %>%
  summarize(Height = mean(value_numeric, na.rm = TRUE)) %>%
  dplyr::select (plant_species = species, Height) %>%
  distinct()

fruit_weight <- Plants %>%
  filter(trait == "FruitWeight") %>%
  group_by(species) %>%
  summarize(MeanFruitWeight = mean(value_numeric, na.rm = TRUE), .groups = 'drop') %>%
  dplyr::select(plant_species = species, MeanFruitWeight) %>%
  distinct()

plants_joined <- Plants %>%
  filter(trait == "fruit_abundance") %>%
  mutate(fruit_abundance = value_numeric) %>%
  left_join(fruit_weight, by = c("species"="plant_species")) 

plants_joined <- plants_joined %>%
  filter(!is.na(MeanFruitWeight)) %>%
  filter(!is.na(fruit_abundance)) %>%
  mutate(CropMass = fruit_abundance * MeanFruitWeight)

crop_mass <- plants_joined %>%
  group_by(species) %>%
  summarize(CropMass = mean(CropMass, na.rm = TRUE), .groups = 'drop') %>%
  dplyr::select(plant_species = species, CropMass)  

traits_plants <- data.frame(plant_species = unique(Plants$species)) %>%
  left_join(fruit_width, by = "plant_species") %>%
  left_join(height, by = "plant_species") %>%
  left_join(crop_mass, by = "plant_species") %>%
  bind_rows(bat_plants)

# write.csv(int_birds, file = here::here("data", "processed", "int_birds.csv"), row.names = F)
# write.csv(int_bat, file = here::here("data", "processed", "int_bats.csv"), row.names = F)
# write.csv(int_nf, file = here::here("data", "processed", "int_nf.csv"), row.names = F)
# 
# write.csv(traits_birds, file = here::here("data", "processed", "traits_birds.csv"), row.names = F)
# write.csv(traits_bat, file = here::here("data", "processed","traits_bats.csv"), row.names = F)
# write.csv(traits_nf, file = here::here("data", "processed", "traits_nf.csv"), row.names = F)
# write.csv(traits_plants, file = here::here("data", "processed", "traits_plants.csv"), row.names = F)

#### Vegetation structure ####

rec <- read.csv("data/raw/Table S4.csv")
con <- read.csv("data/raw/Felicity_data.csv")
vvh <- read.csv("data/raw/verticalVH.csv")

env_variables <- master[c(1:20, 22:62,64),] %>% 
  left_join(select(vvh, Plot_ID, VerticalVH), by = "Plot_ID") %>%
  left_join(select(rec, Plot_ID, Max_tree_height, AGB_all, AGB_wild, sr_all, sr_wild), by = "Plot_ID") %>%
  left_join(select(con, Plot_ID, Forest_1km, Forest_500m, Forest_100m, Distance_forest, Distance_edge, Patch_ha, Cacao_1km, Cacao_Reg1_1km, Cacao_Reg2_1km, Pasture_1km, Pasture_Reg1_1km, Pasture_Reg2_1km), by = "Plot_ID")

env_variables <- env_variables %>%
  mutate(Treatment3 = as.factor(case_when(
    Treatment2 %in% c("active cacao", "active pasture") ~ "regeneration early",
    Treatment2 %in% c("cacao regeneration early", "pasture regeneration early") ~ "regeneration early",
    Treatment2 %in% c("cacao regeneration late", "pasture regeneration late") ~ "regeneration late",
    TRUE ~ as.character(Treatment2) 
  )))

# write.csv(env_variables, file = here::here("data", "processed", "env_variables.csv"), row.names = F)
