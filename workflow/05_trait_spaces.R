##### Trait spaces ####

scores_pol <- readRDS("output/scores_processes.RDS")[[1]]
scores_sd <- readRDS("output/scores_processes.RDS")[[2]]

pca_pol <- readRDS("output/pca_processes.RDS")[[1]]
pca_sd <- readRDS("output/pca_processes.RDS")[[2]]

theme_manuscript <- theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    plot.title = element_text(face = "bold", hjust = 0),
    axis.ticks.length = unit(-0.15, "cm"),
    plot.margin = margin(10, 5, 10, 10)
  )

# Pollination

scores_pol$group <- factor(scores_pol$group, levels = c("Moths", "Bees", "Bat_pol"))

hull_pol <- hulls_groups(scores_pol)

hull_pol$group <- factor(hull_pol$group, levels = c("Moths", "Bees", "Bat_pol"))

centroids_pol <- aggregate(cbind(RC1, RC2) ~ group, data = scores_pol, FUN = mean)

ts_pol <- ggplot(scores_pol, aes(x = RC1, y = RC2, color=group)) +
  geom_point(data = scores_pol, aes(RC1, RC2)) +
  scale_color_manual(values = c('Bat_pol' = "#c994c7",
                                'Bees' =  "#9e9ac8",
                                'Moths' = "#cbc9e2"),
                     labels = c("Bats", "Bees", "Moths")) +
  geom_polygon(data = hull_pol, aes(x = RC1, y = RC2, color = group, fill = group), alpha = 0.4) +
  scale_fill_manual(values = c('Moths' = "#cbc9e2",
                               'Bees' =  "#9e9ac8",
                               'Bat_pol' = "#c994c7"),
                    labels = c("Bats", "Bees", "Moths")) +
  # geom_point(data = centroids_pol, aes(x = RC1, y = RC2, fill = group), 
  #            size = 5, shape = 23, color = "black", stroke = 1) +
  # geom_segment(data = as.data.frame(pca_pol$loadings[,1:2] * 3), aes(x = 0, y = 0, xend = RC1, yend = RC2),
  #              arrow = arrow(), color = "black", linewidth = 1) +
  # geom_text(data = as.data.frame(pca_pol$loadings[, 1:2] * 3), aes(label = c("ProbLength", "CorolLength", "WingSize",  "Height"),
  #                                                                  x = RC1, y = RC2), color = "black", size = 4, vjust = -1, hjust = 0.2, fontface = "bold") +
  # coord_fixed(ratio = 1) +
  theme_manuscript +
  ggtitle("A) Pollination") +
  # theme_classic() +
  # theme(legend.text = element_text(size=16),
  #       legend.position = "none",
  #       axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
  #       # axis.title = element_blank()
  # ) +
  scale_x_continuous(limits = c(floor(range(scores_pol$RC1)[1]), ceiling(range(scores_pol$RC1)[2]))) +
  scale_y_continuous(limits = c(floor(range(scores_pol$RC2)[1]), ceiling(range(scores_pol$RC2)[2]))) +
  labs(x = paste("Trait axis 1 (", round(pca_pol$Vaccounted[2,1]*100, 2), "%)", sep = ""),
       y = paste("Trait axis 2 (", round(pca_pol$Vaccounted[2,2]*100, 2), "%)", sep = ""))
ts_pol

loadings_pol <- ggplot(scores_pol, aes(x = RC1, y = RC2, color=group)) + 
  geom_segment(data = as.data.frame(pca_pol$loadings[,1:2] * 3), aes(x = 0, y = 0, xend = RC1, yend = RC2),
               arrow = arrow(), color = "black", linewidth = 1) +
  geom_text(data = as.data.frame(pca_pol$loadings[, 1:2] * 3), aes(label = c("ProbLength", "CorolLength", "WingSize",  "Height"),
                                                                   x = RC1, y = RC2), color = "black", size = 4, vjust = -1, hjust = 0.2, fontface = "bold") +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits = c(floor(range(scores_pol$RC1)[1]), ceiling(range(scores_pol$RC1)[2]))) +
  scale_y_continuous(limits = c(floor(range(scores_pol$RC2)[1]), ceiling(range(scores_pol$RC2)[2]))) +
  theme_manuscript 

# Seed dispersal

scores_sd$group <- factor(scores_sd$group, levels = c("Bats", "Birds", "NF"))

hull_sd <- hulls_groups(scores_sd)

hull_sd$group <- factor(hull_sd$group, levels = c("Bats", "Birds", "NF"))

centroids_sd <- aggregate(cbind(RC1, RC2) ~ group, data = scores_sd, FUN = mean)

ts_sd <- ggplot(scores_sd, aes(x = RC1, y = RC2, color=group)) +
  geom_point(data = scores_sd, aes(RC1, RC2)) +
  scale_color_manual(values = c('Bats' = "#3182bd",
                                'Birds' =  "#7bccc4",
                                'NF' = "#bdd7e7"),
                     labels = c("Bats", "Birds", "Non-flying/nmammals")) +
  geom_polygon(data = hull_sd, aes(x = RC1, y = RC2, color = group, fill = group), alpha = 0.4) +
  scale_fill_manual(values = c('Bats' = "#3182bd",
                               'Birds' =  "#7bccc4",
                               'NF' = "#bdd7e7"),
                    labels = c("Bats", "Birds", "Non-flying/nmammals")) +
  # geom_point(data = centroids_sd, aes(x = RC1, y = RC2, fill = group), 
  #            size = 5, shape = 23, color = "black", stroke = 1) +
  # geom_segment(data = as.data.frame(pca_sd$loadings[,1:2] * 3), aes(x = 0, y = 0, xend = RC1, yend = RC2),
  #              arrow = arrow(), color = "black", linewidth = 1) +
  # geom_text(data = as.data.frame(pca_sd$loadings[, 1:2] * 3), aes(label = c("BodyMass", "CropMass", "GapeWidth", "FruitWidth", "HandWingIndex", "Height"),
  #                                                                 x = RC1, y = RC2), color = "black", size = 4, vjust = -1, hjust = 0.2, fontface = "bold") +
  # coord_fixed(ratio = 1) +
  theme_manuscript +
  # theme_classic(base_size = 16) +
  # theme(legend.text = element_text(size=16),
  #       axis.title = element_text(size=16),,
  #       plot.title = element_text(face = "bold", hjust = 0),
  #       legend.position = "none",
  #       axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
  #       plot.margin = margin(10, 5, 10, 10),
  #       axis.ticks.length = unit(-0.15, "cm")
  #       # axis.title = element_blank()
  # ) +
  ggtitle("B) Seed dispersal") +
  scale_x_continuous(limits = c(floor(range(scores_sd$RC1)[1]), ceiling(range(scores_sd$RC1)[2]))) +
  scale_y_continuous(limits = c(floor(range(scores_sd$RC2)[1]), ceiling(range(scores_sd$RC2)[2]))) +
  labs(x = paste("Trait axis 1 (", round(pca_sd$Vaccounted[2,1]*100, 2), "%)", sep = ""),
       y = paste("Trait axis 2 (", round(pca_sd$Vaccounted[2,2]*100, 2), "%)", sep = ""))
ts_sd

loadings_sd <- ggplot(scores_sd, aes(x = RC1, y = RC2, color=group)) +
  geom_segment(data = as.data.frame(pca_sd$loadings[,1:2] * 3), aes(x = 0, y = 0, xend = RC1, yend = RC2),
               arrow = arrow(), color = "black", linewidth = 1) +
  geom_text(data = as.data.frame(pca_sd$loadings[, 1:2] * 3), aes(label = c("BodyMass", "CropMass", "GapeWidth", "FruitWidth", "HandWingIndex", "Height"),
                                                                  x = RC1, y = RC2), color = "black", size = 4, vjust = -1, hjust = 0.2, fontface = "bold") +
  coord_fixed(ratio = 1) +
  theme_manuscript +
  scale_x_continuous(limits = c(floor(range(scores_sd$RC1)[1]), ceiling(range(scores_sd$RC1)[2]))) +
  scale_y_continuous(limits = c(floor(range(scores_sd$RC2)[1]), ceiling(range(scores_sd$RC2)[2]))) +
  theme_manuscript

ts <- ts_pol + ts_sd + plot_layout(ncol = 2)
loads <- loadings_pol +  loadings_sd + plot_layout(ncol = 2)

svg("output/Figures/Figure_trait_space.svg", width = 10, height = 5)
ts
dev.off()

svg("output/Figures/Figure_trait_space_loads.svg", width = 10, height = 5)
loads
dev.off()

###separating per recovery stage:
stage_order <- c("regeneration early", "regeneration late", "old-growth forest")
scores_pol$Treatment3 <- factor(scores_pol$Treatment3, levels = stage_order)
scores_sd$Treatment3  <- factor(scores_sd$Treatment3,  levels = stage_order)

stage_labeller <- as_labeller(c(
  "regeneration early" = "Early recovery",
  "regeneration late" = "Late recovery",
  "old-growth forest" = "Old-growth forest"
))

# ---- Pollination ----
scores_pol$group <- factor(scores_pol$group, levels = c("Moths", "Bees", "Bat_pol"))
hull_pol <- hulls_groups_stage(scores_pol)
hull_pol$group <- factor(hull_pol$group, levels = c("Moths", "Bees", "Bat_pol"))

ts_pol <- ggplot(scores_pol, aes(x = RC1, y = RC2, color = group)) +
  geom_point() +
  scale_color_manual(values = c('Bat_pol' = "#c994c7", 'Bees' = "#9e9ac8", 'Moths' = "#cbc9e2"),
                     labels = c("Bats", "Bees", "Moths")) +
  geom_polygon(data = hull_pol, aes(x = RC1, y = RC2, group = group, color = group, fill = group), alpha = 0.4) +
  scale_fill_manual(values = c('Moths' = "#cbc9e2", 'Bees' = "#9e9ac8", 'Bat_pol' = "#c994c7"),
                    labels = c("Bats", "Bees", "Moths")) +
  facet_wrap(~Treatment3, ncol = 3, labeller = stage_labeller) +
  theme_manuscript +
  theme(legend.position = "bottom") +   # you'll likely want a legend now, since color = group across panels
  ggtitle("A) Pollination") +
  scale_x_continuous(limits = c(floor(range(scores_pol$RC1)[1]), ceiling(range(scores_pol$RC1)[2]))) +
  scale_y_continuous(limits = c(floor(range(scores_pol$RC2)[1]), ceiling(range(scores_pol$RC2)[2]))) +
  labs(x = paste("Trait axis 1 (", round(pca_pol$Vaccounted[2,1]*100, 2), "%)", sep = ""),
       y = paste("Trait axis 2 (", round(pca_pol$Vaccounted[2,2]*100, 2), "%)", sep = ""))

# ---- Seed dispersal ----
scores_sd$group <- factor(scores_sd$group, levels = c("Bats", "Birds", "NF"))
hull_sd <- hulls_groups_stage(scores_sd)
hull_sd$group <- factor(hull_sd$group, levels = c("Bats", "Birds", "NF"))

ts_sd <- ggplot(scores_sd, aes(x = RC1, y = RC2, color = group)) +
  geom_point() +
  scale_color_manual(values = c('Bats' = "#3182bd", 'Birds' = "#7bccc4", 'NF' = "#bdd7e7"),
                     labels = c("Bats", "Birds", "Non-flying mammals")) +
  geom_polygon(data = hull_sd, aes(x = RC1, y = RC2, group = group, color = group, fill = group), alpha = 0.4) +
  scale_fill_manual(values = c('Bats' = "#3182bd", 'Birds' = "#7bccc4", 'NF' = "#bdd7e7"),
                    labels = c("Bats", "Birds", "Non-flying mammals")) +
  facet_wrap(~Treatment3, ncol = 3, labeller = stage_labeller) +
  theme_manuscript +
  theme(legend.position = "bottom") +
  ggtitle("B) Seed dispersal") +
  scale_x_continuous(limits = c(floor(range(scores_sd$RC1)[1]), ceiling(range(scores_sd$RC1)[2]))) +
  scale_y_continuous(limits = c(floor(range(scores_sd$RC2)[1]), ceiling(range(scores_sd$RC2)[2]))) +
  labs(x = paste("Trait axis 1 (", round(pca_sd$Vaccounted[2,1]*100, 2), "%)", sep = ""),
       y = paste("Trait axis 2 (", round(pca_sd$Vaccounted[2,2]*100, 2), "%)", sep = ""))

ts <- ts_pol / ts_sd  # stacked (2 rows x 3 cols each), or use ts_pol + ts_sd for side by side with plot_layout(ncol=1)

svg("output/Figures/Figure_trait_spaces_per_stage.svg", width = 13, height = 9)
ts
dev.off()
