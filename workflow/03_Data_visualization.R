library(ggplot2)

#### Functional trait spaces ####

## Bats

hull_bats <- hulls_hbt(scores_bats)

ts_bats <- ggplot(scores_bats, aes(x = RC1, y = RC2, color=Treatment3)) +
  geom_polygon(data = hull_bats, aes(x = RC1, y = RC2, fill = Treatment3)) +
  scale_fill_manual(values = c('regeneration early' = rgb(0.89,0.15,0.21,0.3),
                               'regeneration late' =  rgb(1,0.75,0,0.3),
                               'old-growth forest' = rgb(0,0.4,0,0.3)),
                    labels = c("Old-growth forest", "Late regeneration", "Early regeneration")) +
  geom_point(data = scores_bats, aes(RC1, RC2)) +
  scale_color_manual(values = c('regeneration early' = rgb(0.89,0.15,0.21,0.5),
                                'regeneration late' =  rgb(1,0.75,0,0.5),
                                'old-growth forest' = rgb(0,0.4,0,0.5)),
                     labels = c("Old-growth forest", "Late regeneration", "Early regeneration")) +
  geom_point(aes(mean(scores_bats[Treatment3=="regeneration early", "RC1"]), mean(scores_bats[Treatment3=="regeneration early", "RC2"])), size = 5, shape = 23, fill =rgb(0.89,0.15,0.21,0.5), color = rgb(0.89,0.15,0.21,0.5)) +
  geom_point(aes(mean(scores_bats[Treatment3=="regeneration late", "RC1"]), mean(scores_bats[Treatment3=="regeneration late", "RC2"])), size = 5, shape = 23, fill = rgb(1, 0.75,0, 0.5), color = rgb(1, 0.75,0, 0.5)) +
  geom_point(aes(mean(scores_bats[Treatment3=="old-growth forest", "RC1"]), mean(scores_bats[Treatment3=="old-growth forest", "RC2"])), size = 5, shape = 23, fill = rgb(0,0.4,0,0.5), color = rgb(0,0.4,0,0.5)) +
  geom_segment(data = as.data.frame(pca_bats$loadings[,1:2] * 3), aes(x = 0, y = 0, xend = RC1, yend = RC2),
               arrow = arrow(), color = "black", linewidth = 1) +
  geom_text(data = as.data.frame(pca_bats$loadings[, 1:2] * 3), aes(label = c("BodyMass", "CropMass", "GapeWidth", "FruitWidth", "HandWingIndex", "Height"),
                                                                   x = RC1, y = RC2), color = "black", size = 4, vjust = -1, hjust = 0.2, fontface = "bold") +
  # geom_text(data =scores_bats[scores_bats$Treatment3 != "old-growth forest" & scores_bats$interaction %in% plant_animal_network$interaction,],  aes(label = interaction), check_overlap = FALSE, size = 2.5, alpha = 1, color = "black")  +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme(legend.text = element_text(size=16),
        legend.position = "none",
        axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        # axis.title = element_blank()
        ) +
  scale_x_continuous(limits = c(-2.5, 4), breaks = c(-2,-1,0,1,2,3,4)) +
  scale_y_continuous(limits = c(-2.5, 4), breaks = c(-2,-1,0,1,2,3,4)) +
labs(x = paste("Trait axis 1 (", round(pca_bats$Vaccounted[2,1]*100, 2), "%)", sep = ""),
     y = paste("Trait axis 2 (", round(pca_bats$Vaccounted[2,2]*100, 2), "%)", sep = ""))
ts_bats

## Birds

hull_birds <- hulls_hbt(scores_birds)

ts_birds <- ggplot(scores_birds, aes(x = RC1, y = RC2, color=Treatment3)) +
  geom_polygon(data = hull_birds, aes(x = RC1, y = RC2, fill = Treatment3)) +
  scale_fill_manual(values = c('regeneration early' = rgb(0.89,0.15,0.21,0.3),
                               'regeneration late' =  rgb(1,0.75,0,0.3),
                               'old-growth forest' = rgb(0,0.4,0,0.3)),
                    labels = c("Old-growth forest", "Late regeneration", "Early regeneration")) +
  geom_point(data = scores_birds, aes(RC1, RC2)) +
  scale_color_manual(values = c('regeneration early' = rgb(0.89,0.15,0.21,0.5),
                                'regeneration late' =  rgb(1,0.75,0,0.5),
                                'old-growth forest' = rgb(0,0.4,0,0.5)),
                     labels = c("Old-growth forest", "Late regeneration", "Early regeneration")) +
  geom_point(aes(mean(scores_birds[Treatment3=="regeneration early", "RC1"]), mean(scores_birds[Treatment3=="regeneration early", "RC2"])), size = 5, shape = 23, fill =rgb(0.89,0.15,0.21,0.5), color = rgb(0.89,0.15,0.21,0.5)) +
  geom_point(aes(mean(scores_birds[Treatment3=="regeneration late", "RC1"]), mean(scores_birds[Treatment3=="regeneration late", "RC2"])), size = 5, shape = 23, fill = rgb(1, 0.75,0, 0.5), color = rgb(1, 0.75,0, 0.5)) +
  geom_point(aes(mean(scores_birds[Treatment3=="old-growth forest", "RC1"]), mean(scores_birds[Treatment3=="old-growth forest", "RC2"])), size = 5, shape = 23, fill = rgb(0,0.4,0,0.5), color = rgb(0,0.4,0,0.5)) +
  geom_segment(data = as.data.frame(pca_birds$loadings[,1:2] * 3), aes(x = 0, y = 0, xend = RC1, yend = RC2),
               arrow = arrow(), color = "black", linewidth = 1) +
  geom_text(data = as.data.frame(pca_birds$loadings[, 1:2] * 3), aes(label = c("BodyMass", "CropMass", "BeakWidth", "FruitWidth", "HandWingIndex", "Height"),
                                                                    x = RC1, y = RC2), color = "black", size = 4, vjust = -1, hjust = 0.2, fontface = "bold") +
  # geom_text(data =scores_birds[scores_birds$Treatment3 != "old-growth forest" & scores_birds$interaction %in% plant_animal_network$interaction,],  aes(label = interaction), check_overlap = FALSE, size = 2.5, alpha = 1, color = "black")  +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme(legend.text = element_text(size=16),
        legend.position = "none",
        axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        # axis.title = element_blank()
        ) +
  scale_x_continuous(limits = c(-2.5, 4), breaks = c(-2,-1,0,1,2,3,4)) +
  scale_y_continuous(limits = c(-2.5, 4), breaks = c(-2,-1,0,1,2,3,4)) +
  labs(x = paste("Trait axis 1 (", round(pca_birds$Vaccounted[2,1]*100, 2), "%)", sep = ""),
       y = paste("Trait axis 2 (", round(pca_birds$Vaccounted[2,2]*100, 2), "%)", sep = ""))
ts_birds

## Non-flying mammals

hull_nf <- hulls_hbt(scores_nf)

ts_nf <- ggplot(scores_nf, aes(x = RC1, y = RC2, color=Treatment3)) +
  geom_polygon(data = hull_nf, aes(x = RC1, y = RC2, fill = Treatment3)) +
  scale_fill_manual(values = c('regeneration early' = rgb(0.89,0.15,0.21,0.3),
                               'regeneration late' =  rgb(1,0.75,0,0.3),
                               'old-growth forest' = rgb(0,0.4,0,0.3)),
                    labels = c("Old-growth forest", "Late regeneration", "Early regeneration")) +
  geom_point(data = scores_nf, aes(RC1, RC2)) +
  scale_color_manual(values = c('regeneration early' = rgb(0.89,0.15,0.21,0.5),
                                'regeneration late' =  rgb(1,0.75,0,0.5),
                                'old-growth forest' = rgb(0,0.4,0,0.5)),
                     labels = c("Old-growth forest", "Late regeneration", "Early regeneration")) +
  geom_point(aes(mean(scores_nf[Treatment3=="regeneration early", "RC1"]), mean(scores_nf[Treatment3=="regeneration early", "RC2"])), size = 5, shape = 23, fill =rgb(0.89,0.15,0.21,0.5), color = rgb(0.89,0.15,0.21,0.5)) +
  geom_point(aes(mean(scores_nf[Treatment3=="regeneration late", "RC1"]), mean(scores_nf[Treatment3=="regeneration late", "RC2"])), size = 5, shape = 23, fill = rgb(1, 0.75,0, 0.5), color = rgb(1, 0.75,0, 0.5)) +
  geom_point(aes(mean(scores_nf[Treatment3=="old-growth forest", "RC1"]), mean(scores_nf[Treatment3=="old-growth forest", "RC2"])), size = 5, shape = 23, fill = rgb(0,0.4,0,0.5), color = rgb(0,0.4,0,0.5)) +
  geom_segment(data = as.data.frame(pca_nf$loadings[,1:2] * 3), aes(x = 0, y = 0, xend = RC1, yend = RC2),
               arrow = arrow(), color = "black", linewidth = 1) +
  geom_text(data = as.data.frame(pca_nf$loadings[, 1:2] * 3), aes(label = c("BodyMass", "CropMass", "GapeWidth", "FruitWidth"),
                                                                    x = RC1, y = RC2), color = "black", size = 4, vjust = -1, hjust = 0.2, fontface = "bold") +
  # geom_text(data =scores_nf[scores_nf$Treatment3 != "old-growth forest" & scores_nf$interaction %in% plant_animal_network$interaction,],  aes(label = interaction), check_overlap = FALSE, size = 2.5, alpha = 1, color = "black")  +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme(legend.text = element_text(size=16),
        legend.position = "none",
        axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        # axis.title = element_blank()
        ) +
  scale_x_continuous(limits = c(-2.5, 4), breaks = c(-2,-1,0,1,2,3,4)) +
  scale_y_continuous(limits = c(-2.5, 4), breaks = c(-2,-1,0,1,2,3,4)) +
  labs(x = paste("Trait axis 1 (", round(pca_nf$Vaccounted[2,1]*100, 2), "%)", sep = ""),
       y = paste("Trait axis 2 (", round(pca_nf$Vaccounted[2,2]*100, 2), "%)", sep = ""))
ts_nf


#### Scatter plots

plot(FDBats~RegTime, model_df)

m_bats <- lm(FDBats ~ Treatment3,
         data = model_df)
summary(m_bats)

plot(FDBirds~RegTime, model_df)

m_birds <- lm(FDBirds ~ Treatment3,
             data = model_df)
summary(m_birds)

plot(FDNf~RegTime, model_df)

m_nf <- lm(FDNf ~ Treatment3,
             data = model_df)
summary(m_nf)
