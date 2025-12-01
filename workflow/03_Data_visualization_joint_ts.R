library(ggplot2)

#### Functional trait spaces ####

## Bats

hull_bats <- hulls_hbt(scores_sd[scores_sd$group == "Bats",])

centroids_bats <- scores_sd %>%
  filter(group == "Bats") %>%
  group_by(Treatment3) %>%
  summarise(RC1 = mean(RC1, na.rm = TRUE),
            RC2 = mean(RC2, na.rm = TRUE), .groups = "drop") 

ts_bats <- ggplot(scores_sd, aes(x = RC1, y = RC2, color=Treatment3)) +
  geom_polygon(data = hull_bats, aes(x = RC1, y = RC2, fill = Treatment3)) +
  scale_fill_manual(values = c('regeneration early' = rgb(0.89,0.15,0.21,0.3),
                               'regeneration late' =  rgb(1,0.75,0,0.3),
                               'old-growth forest' = rgb(0,0.4,0,0.3)),
                    labels = c("Old-growth forest", "Late regeneration", "Early regeneration")) +
  geom_point(data = scores_sd[scores_sd$group == "Bats",], aes(RC1, RC2)) +
  scale_color_manual(values = c('regeneration early' = rgb(0.89,0.15,0.21,0.5),
                                'regeneration late' =  rgb(1,0.75,0,0.5),
                                'old-growth forest' = rgb(0,0.4,0,0.5)),
                     labels = c("Old-growth forest", "Late regeneration", "Early regeneration")) +
  geom_point(data = centroids_bats, aes(RC1, RC2, fill = Treatment3, color = Treatment3), size = 5, shape = 23) +  
  # geom_segment(data = as.data.frame(pca_sd$loadings[,1:2] * 3), aes(x = 0, y = 0, xend = RC1, yend = RC2),
  #              arrow = arrow(), color = "black", linewidth = 1) +
  # geom_text(data = as.data.frame(pca_sd$loadings[, 1:2] * 3), aes(label = c("BodyMass", "CropMass", "GapeWidth", "FruitWidth", "HandWingIndex", "Height"),
  #                                                                 x = RC1, y = RC2), color = "black", size = 4, vjust = -1, hjust = 0.2, fontface = "bold") +
  geom_text(data = scores_sd[scores_sd$group == "Bats",],  aes(label = interaction), check_overlap = FALSE, size = 2.5, alpha = 1, color = "black")  +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme(legend.text = element_text(size=16),
        legend.position = "none",
        axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        # axis.title = element_blank()
  ) +
  scale_x_continuous(limits = c(-3.1, 3.5), breaks = c(-2,0,2)) +
  scale_y_continuous(limits = c(-2.3, 4), breaks = c(-2,0,2,4)) +
  labs(x = paste("Trait axis 1 (", round(pca_sd$Vaccounted[2,1]*100, 2), "%)", sep = ""),
       y = paste("Trait axis 2 (", round(pca_sd$Vaccounted[2,2]*100, 2), "%)", sep = ""))
ts_bats

## Birds

hull_birds <- hulls_hbt(scores_sd[scores_sd$group == "Birds",])

ts_birds <- ggplot(scores_sd, aes(x = RC1, y = RC2, color=Treatment3)) +
  geom_polygon(data = hull_birds, aes(x = RC1, y = RC2, fill = Treatment3)) +
  scale_fill_manual(values = c('regeneration early' = rgb(0.89,0.15,0.21,0.3),
                               'regeneration late' =  rgb(1,0.75,0,0.3),
                               'old-growth forest' = rgb(0,0.4,0,0.3)),
                    labels = c("Old-growth forest", "Late regeneration", "Early regeneration")) +
  geom_point(data = scores_sd[scores_sd$group == "Birds",], aes(RC1, RC2)) +
  scale_color_manual(values = c('regeneration early' = rgb(0.89,0.15,0.21,0.5),
                                'regeneration late' =  rgb(1,0.75,0,0.5),
                                'old-growth forest' = rgb(0,0.4,0,0.5)),
                     labels = c("Old-growth forest", "Late regeneration", "Early regeneration")) +
  geom_point(aes(mean(scores_sd[Treatment3=="regeneration early" & group == "Birds", "RC1"]), mean(scores_sd[Treatment3=="regeneration early" & group == "Birds", "RC2"])), size = 5, shape = 23, fill =rgb(0.89,0.15,0.21,0.5), color = rgb(0.89,0.15,0.21,0.5)) +
  geom_point(aes(mean(scores_sd[Treatment3=="regeneration late" & group == "Birds", "RC1"]), mean(scores_sd[Treatment3=="regeneration late" & group == "Birds", "RC2"])), size = 5, shape = 23, fill = rgb(1, 0.75,0, 0.5), color = rgb(1, 0.75,0, 0.5)) +
  geom_point(aes(mean(scores_sd[Treatment3=="old-growth forest" & group == "Birds", "RC1"]), mean(scores_sd[Treatment3=="old-growth forest" & group == "Birds", "RC2"])), size = 5, shape = 23, fill = rgb(0,0.4,0,0.5), color = rgb(0,0.4,0,0.5)) +
  # geom_segment(data = as.data.frame(pca_sd$loadings[,1:2] * 3), aes(x = 0, y = 0, xend = RC1, yend = RC2),
  #              arrow = arrow(), color = "black", linewidth = 1) +
  # geom_text(data = as.data.frame(pca_sd$loadings[, 1:2] * 3), aes(label = c("BodyMass", "CropMass", "GapeWidth", "FruitWidth", "HandWingIndex", "Height"),
  #                                                                 x = RC1, y = RC2), color = "black", size = 4, vjust = -1, hjust = 0.2, fontface = "bold") +
  geom_text(data = scores_sd[scores_sd$group == "Birds",],  aes(label = interaction), check_overlap = FALSE, size = 2.5, alpha = 1, color = "black")  +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme(legend.text = element_text(size=16),
        legend.position = "none",
        axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        # axis.title = element_blank()
  ) +
  scale_x_continuous(limits = c(-3.1, 3.5), breaks = c(-2,0,2)) +
  scale_y_continuous(limits = c(-2.3, 4), breaks = c(-2,0,2,4)) +
  labs(x = paste("Trait axis 1 (", round(pca_sd$Vaccounted[2,1]*100, 2), "%)", sep = ""),
       y = paste("Trait axis 2 (", round(pca_sd$Vaccounted[2,2]*100, 2), "%)", sep = ""))
ts_birds

## Non-flying mammals

hull_nf <- hulls_hbt(scores_sd[scores_sd$group == "NF",])

ts_nf <- ggplot(scores_sd, aes(x = RC1, y = RC2, color=Treatment3)) +
  geom_polygon(data = hull_nf, aes(x = RC1, y = RC2, fill = Treatment3)) +
  scale_fill_manual(values = c('regeneration early' = rgb(0.89,0.15,0.21,0.3),
                               'regeneration late' =  rgb(1,0.75,0,0.3),
                               'old-growth forest' = rgb(0,0.4,0,0.3)),
                    labels = c("Old-growth forest", "Late regeneration", "Early regeneration")) +
  geom_point(data = scores_sd[scores_sd$group == "NF",], aes(RC1, RC2)) +
  scale_color_manual(values = c('regeneration early' = rgb(0.89,0.15,0.21,0.5),
                                'regeneration late' =  rgb(1,0.75,0,0.5),
                                'old-growth forest' = rgb(0,0.4,0,0.5)),
                     labels = c("Old-growth forest", "Late regeneration", "Early regeneration")) +
  geom_point(aes(mean(scores_sd[Treatment3=="regeneration early" & group == "NF", "RC1"]), mean(scores_sd[Treatment3=="regeneration early" & group == "NF", "RC2"])), size = 5, shape = 23, fill =rgb(0.89,0.15,0.21,0.5), color = rgb(0.89,0.15,0.21,0.5)) +
  geom_point(aes(mean(scores_sd[Treatment3=="regeneration late" & group == "NF", "RC1"]), mean(scores_sd[Treatment3=="regeneration late" & group == "NF", "RC2"])), size = 5, shape = 23, fill = rgb(1, 0.75,0, 0.5), color = rgb(1, 0.75,0, 0.5)) +
  geom_point(aes(mean(scores_sd[Treatment3=="old-growth forest" & group == "NF", "RC1"]), mean(scores_sd[Treatment3=="old-growth forest" & group == "NF", "RC2"])), size = 5, shape = 23, fill = rgb(0,0.4,0,0.5), color = rgb(0,0.4,0,0.5)) +
  # geom_segment(data = as.data.frame(pca_sd$loadings[,1:2] * 3), aes(x = 0, y = 0, xend = RC1, yend = RC2),
  #              arrow = arrow(), color = "black", linewidth = 1) +
  # geom_text(data = as.data.frame(pca_sd$loadings[, 1:2] * 3), aes(label = c("BodyMass", "CropMass", "GapeWidth", "FruitWidth", "HandWingIndex", "Height"),
  #                                                                 x = RC1, y = RC2), color = "black", size = 4, vjust = -1, hjust = 0.2, fontface = "bold") +
  geom_text(data = scores_sd[scores_sd$group == "NF",],  aes(label = interaction), check_overlap = FALSE, size = 2.5, alpha = 1, color = "black")  +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme(legend.text = element_text(size=16),
        legend.position = "none",
        axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        # axis.title = element_blank()
  ) +
  scale_x_continuous(limits = c(-3.1, 3.5), breaks = c(-2,0,2)) +
  scale_y_continuous(limits = c(-2.3, 4), breaks = c(-2,0,2,4)) +
  labs(x = paste("Trait axis 1 (", round(pca_sd$Vaccounted[2,1]*100, 2), "%)", sep = ""),
       y = paste("Trait axis 2 (", round(pca_sd$Vaccounted[2,2]*100, 2), "%)", sep = ""))
ts_nf

## ALL SD

hull_sd <- hulls_hbt(scores_sd)

ts_sd <- ggplot(scores_sd, aes(x = RC1, y = RC2, color=Treatment3)) +
  geom_polygon(data = hull_sd, aes(x = RC1, y = RC2, fill = Treatment3)) +
  scale_fill_manual(values = c('regeneration early' = rgb(0.89,0.15,0.21,0.3),
                               'regeneration late' =  rgb(1,0.75,0,0.3),
                               'old-growth forest' = rgb(0,0.4,0,0.3)),
                    labels = c("Old-growth forest", "Late regeneration", "Early regeneration")) +
  # geom_point(data = scores_sd, aes(RC1, RC2)) +
  scale_color_manual(values = c('regeneration early' = rgb(0.89,0.15,0.21,0.5),
                                'regeneration late' =  rgb(1,0.75,0,0.5),
                                'old-growth forest' = rgb(0,0.4,0,0.5)),
                     labels = c("Old-growth forest", "Late regeneration", "Early regeneration")) +
 
  # geom_point(aes(mean(scores_sd[Treatment3=="regeneration early", "RC1"]), mean(scores_nf[Treatment3=="regeneration early", "RC2"])), size = 5, shape = 23, fill =rgb(0.89,0.15,0.21,0.5), color = rgb(0.89,0.15,0.21,0.5)) +
  # geom_point(aes(mean(scores_sd[Treatment3=="regeneration late", "RC1"]), mean(scores_nf[Treatment3=="regeneration late", "RC2"])), size = 5, shape = 23, fill = rgb(1, 0.75,0, 0.5), color = rgb(1, 0.75,0, 0.5)) +
  # geom_point(aes(mean(scores_sd[Treatment3=="old-growth forest", "RC1"]), mean(scores_nf[Treatment3=="old-growth forest", "RC2"])), size = 5, shape = 23, fill = rgb(0,0.4,0,0.5), color = rgb(0,0.4,0,0.5)) +
  geom_segment(data = as.data.frame(pca_sd$loadings[,1:2] * 3), aes(x = 0, y = 0, xend = RC1, yend = RC2),
               arrow = arrow(), color = "black", linewidth = 1) +
  geom_text(data = as.data.frame(pca_sd$loadings[, 1:2] * 3), aes(label = c("BodyMass", "CropMass", "GapeWidth", "FruitWidth", "HandWingIndex", "Height"),
                                                                     x = RC1, y = RC2), color = "black", size = 4, vjust = -1, hjust = 0.2, fontface = "bold") +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme(legend.text = element_text(size=16),
        legend.position = "none",
        axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        # axis.title = element_blank()
  ) +
  scale_x_continuous(limits = c(-3.1, 3.5), breaks = c(-2,0,2)) +
  scale_y_continuous(limits = c(-2.3, 4), breaks = c(-2,0,2,4)) +
  labs(x = paste("Trait axis 1 (", round(pca_sd$Vaccounted[2,1]*100, 2), "%)", sep = ""),
       y = paste("Trait axis 2 (", round(pca_sd$Vaccounted[2,2]*100, 2), "%)", sep = ""))
ts_sd

#### Scatter plots

# FD

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

# FC1

plot(FC1Bats~RegTime, model_df)

m_bats <- lm(FC1Bats ~ Treatment3,
             data = model_df)
summary(m_bats)

plot(FC1Birds~RegTime, model_df)

m_birds <- lm(FC1Birds ~ Treatment3,
              data = model_df)
summary(m_birds)

plot(FC1Nf~RegTime, model_df)

m_nf <- lm(FC1Nf ~ Treatment3,
           data = model_df)
summary(m_nf)

# FC2

plot(FC2Bats~RegTime, model_df)

m_bats <- lm(FC2Bats ~ Treatment3,
             data = model_df)
summary(m_bats)

plot(FC2Birds~RegTime, model_df)

m_birds <- lm(FC2Birds ~ Treatment3,
              data = model_df)
summary(m_birds)

plot(FC2Nf~RegTime, model_df)

m_nf <- lm(FC2Nf ~ Treatment3,
           data = model_df)
summary(m_nf)
