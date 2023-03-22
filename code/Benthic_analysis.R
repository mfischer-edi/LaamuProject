# Laamu Coral Reef monitoring project
# Script for analysis of spatiotemporal changes in the benthic community composition,
# and analysis of anthropogenic and non-human reef impacts.

# Written by Mara Fischer
# E-mail: mf555@exeter.ac.uk
# 02/09/2022

# Packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(DHARMa)
library(lme4)
library(lmerTest)
library(ggeffects)
library(MuMIn)
library(vegan)

# Create theme for plots
theme_diss <- function(){            
  theme_bw()+                          
    theme(text = element_text(family = "sans"),
          axis.text = element_text(size = 16), 
          axis.title = element_text(size = 20),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 20, vjust = 1, hjust = 0),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16, face = "bold"),
          legend.position = "bottom",
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent",
                                           size = 2, linetype = "blank"))
}

# Load data ----

# 2019 and 2021/22 benthic data for temporal analysis
temp_an <- read.csv("data/all_perc_data_final.csv", header = T) # full data
temp_an_dep <- read.csv("data/perc_site_dep.csv", header = T)   # data aggr. by depth
perc_site <- read.csv("data/perc_site.csv", header = T)         # data aggr. by site
   
# All 2021/22 benthic data 
all22sites <- read.csv("data/all22sites.csv", header = T)  # full data
all22sites_aggr <- read.csv("data/all22sites_aggr.csv", header = T) # (aggregated categories)

# Coral cover by site
coral22 <- all22sites_aggr %>% 
  group_by(Site.name) %>% 
  summarise(Hard.coral = mean(Hard.coral)) %>% 
  ungroup()

write.csv(coral22, "data/coral22_site.csv")

# % change over time in coral cover and impacts
change <- read.csv("data/coral_rate_of_change.csv", header = T)

# Benthic data aggr. by site with impacts by type
all_data <- read.csv("data/all_data_combined.csv", header = T)  

# Raw impacts data
impacts22 <- read.csv("data/impacts_data2.csv", header = T) # full 2022 raw impacts data
impacts19 <- read.csv("data/Impacts19_raw.csv", header = T) # full 2019 raw impacts data

# ANALYSIS ----

# Benthic community ----

## Year NMDS ----

# Use data agrgegated by year and site

# Subset dataframe, on which to base the ordination:
data1 <- perc_site[,6:16]
# Descriptive columns
data2 <- perc_site[,2:3]

set.seed(1506)
# Because the final result depends on the initial random placement of the points, 
# I`ll set a seed to make the results reproducible (get same outcome every time)

nmds2 <- metaMDS(data1, k = 2, trymax = 100, trace = F, autotransform = FALSE, distance="bray")
nmds2

stressplot(nmds2)
plot(nmds2)

ordiplot(nmds2, type = "n")
orditorp(nmds2, display = "species", col = "red", air = 0.01)
orditorp(nmds2, display = "site", cex = 1.1, air = 0.01)

data.scores2 <- as.data.frame(nmds2$points)
nmds_all2 <- cbind(data2, data.scores2)

species.scores <- as.data.frame(scores(nmds2, "species"))
species.scores$category <- rownames(species.scores)

# Make hulls
grp.a <- nmds_all2[nmds_all2$year == "2019", ][chull(nmds_all2[nmds_all2$year == "2019", c("MDS1", "MDS2")]), ]  
grp.b <- nmds_all2[nmds_all2$year == "2021/22", ][chull(nmds_all2[nmds_all2$year == "2021/22", c("MDS1", "MDS2")]), ]  
hull.data <- rbind(grp.a, grp.b)  #combine grp.a and grp.b
hull.data

# Plot with ggplot
(nmds_plot <- ggplot() + 
    geom_polygon(data = hull.data, aes(x = MDS1, y = MDS2, fill = year, group = year), alpha = 0.30) + # add the convex hulls
    geom_text(data = species.scores, aes(x = NMDS1, y = NMDS2, label = category), 
              alpha=0.5, size = 4.5) +  # add the species labels
    geom_point(data = nmds_all2, aes(x = MDS1, y = MDS2, shape = year, colour = year),size = 3) + # add the point markers
    #geom_text(data = nmds_all2, aes(x = NMDS1, y = NMDS2, label = survey.id), size = 1, vjust = 0) +  # add the site labels
    #annotate("text", x = -0.68, y = -0.76, label = "Ecklonia radiata", size = 3) + 
    labs(x = "NMDS1", y = "NDMS2") +
    scale_x_continuous(limits = c(-1.2, 0.8)) +
    scale_colour_manual(values = c("2019" = "#ffba49",
                                   "2021/22" = "#20a39e")) +
    scale_fill_manual(values = c("2019" = "#ffba49",
                                 "2021/22" = "#20a39e"))+
    coord_equal() +
    theme_diss() +
    theme(panel.border = element_rect(fill = NA),
          legend.title = element_blank(),
          legend.position = c(0.85,0.9)))

ggsave(nmds_plot, filename = "outputs/nmds_year.png", width = 10, height = 8)

## Year ANOSIM & SIMPER ----

group <- substr(data2$year, 1,38)
group[group=="1"] <- "2019"
group[group=="2"] <- "2021/22"

group.fac <- factor(group, levels = c("2019", "2021/22"))

# setting the seed means the 'random' outcome will be reproducible.
set.seed(123)
anosim(data1, data2$year, permutations = 999, distance = "bray", strata = NULL,
       parallel = getOption("mc.cores"))

group = rep(c("2019", "2021/22"), 38)
group.fac <- factor(group, levels = c("2019", "2021/22"))

# Year SIMPER

simper(data1, data2$year, permutations = 999, trace = FALSE)

## Depth NMDS ----

# Use all 2022 data but same categories

# Reorder depths
all22sites_aggr$Depth <- factor(all22sites_aggr$Depth,
                                levels = c("5m", "10m", "15m"),
                                labels = c("5m", "10m", "15m"))


# Subset dataframe, on which to base the ordination:
data1 <- all22sites_aggr[,6:16]
# Descriptive columns
data2 <- all22sites_aggr[,2:3]

set.seed(1506)
# Because the final result depends on the initial random placement of the points, 
# I`ll set a seed to make the results reproducible (get same outcome every time)

nmds_dep <- metaMDS(data1, k = 2, trymax = 100, trace = F, autotransform = FALSE, distance="bray")
nmds_dep

stressplot(nmds_dep)
plot(nmds_dep)

ordiplot(nmds_dep, type = "n")
orditorp(nmds_dep, display = "species", col = "red", air = 0.01)
orditorp(nmds_dep, display = "site", cex = 1.1, air = 0.01)

data.scores2 <- as.data.frame(nmds_dep$points)
nmds_all2 <- cbind(data2, data.scores2)

species.scores <- as.data.frame(scores(nmds_dep, "species"))
species.scores$category <- rownames(species.scores)

# Make hulls
grp.a <- nmds_all2[nmds_all2$Depth == "5m", ][chull(nmds_all2[nmds_all2$Depth == "5m", c("MDS1", "MDS2")]), ]  
grp.b <- nmds_all2[nmds_all2$Depth == "10m", ][chull(nmds_all2[nmds_all2$Depth== "10m", c("MDS1", "MDS2")]), ]  
grp.c <- nmds_all2[nmds_all2$Depth == "15m", ][chull(nmds_all2[nmds_all2$Depth== "15m", c("MDS1", "MDS2")]), ]  

hull.data <- rbind(grp.a, grp.b, grp.c)  #combine
hull.data

# Plot with ggplot
(nmds_plot_dep <- ggplot() + 
    geom_point(data = nmds_all2, aes(x = MDS1, y = MDS2, shape = Depth, colour = Depth),size = 3) + # add the point markers
    geom_polygon(data = hull.data, aes(x = MDS1, y = MDS2, fill = Depth, group = Depth), alpha = 0.40) + # add the convex hulls
    geom_text(data = species.scores, aes(x = NMDS1, y = NMDS2, label = category), 
              alpha=0.5, size = 4.5) +  # add the species labels
    labs(x = "NMDS1", y = "NDMS2", colour = "Depth", shape = "Depth") +
    scale_colour_manual(values = c('#2c7fb8','#7fcdbb','#edf8b1')) +
    scale_fill_manual(values = c('#2c7fb8','#7fcdbb','#edf8b1'))+
    coord_equal() +
    theme_diss() +
    theme(panel.border = element_rect(fill = NA),
          legend.title = element_blank(),
          legend.position = c(0.85,0.9)))

ggsave(nmds_plot_dep, filename = "outputs/nmds_depth.png", width = 10, height = 8)

## Reef type NMDS ----

# Use all 2022 data but same categories

# Subset dataframe, on which to base the ordination:
data1 <- all22sites_aggr[,6:16]
# Descriptive columns
data2 <- all22sites_aggr[,4:5]

set.seed(1506)
# Because the final result depends on the initial random placement of the points, 
# I`ll set a seed to make the results reproducible (get same outcome every time)

nmds_rt <- metaMDS(data1, k = 2, trymax = 100, trace = F, autotransform = FALSE, distance="bray")
nmds_rt

stressplot(nmds_rt)
plot(nmds_rt)

ordiplot(nmds_rt, type = "n")
orditorp(nmds_rt, display = "species", col = "red", air = 0.01)
orditorp(nmds_rt, display = "site", cex = 1.1, air = 0.01)

data.scores2 <- as.data.frame(nmds_rt$points)
nmds_all2 <- cbind(data2, data.scores2)

species.scores <- as.data.frame(scores(nmds_rt, "species"))
species.scores$category <- rownames(species.scores)

# Make hulls
grp.a <- nmds_all2[nmds_all2$Reef.type == "Channel", ][chull(nmds_all2[nmds_all2$Reef.type == "Channel", c("MDS1", "MDS2")]), ]  
grp.b <- nmds_all2[nmds_all2$Reef.type == "Outer", ][chull(nmds_all2[nmds_all2$Reef.type == "Outer", c("MDS1", "MDS2")]), ]  
grp.c <- nmds_all2[nmds_all2$Reef.type == "Inner", ][chull(nmds_all2[nmds_all2$Reef.type == "Inner", c("MDS1", "MDS2")]), ]  
grp.d <- nmds_all2[nmds_all2$Reef.type == "Thila", ][chull(nmds_all2[nmds_all2$Reef.type == "Thila", c("MDS1", "MDS2")]), ]
hull.data <- rbind(grp.a, grp.b, grp.c, grp.d)  #combine
hull.data

# Plot with ggplot
(nmds_plot_rt <- ggplot() + 
    geom_point(data = nmds_all2, aes(x = MDS1, y = MDS2, shape = Reef.type, colour = Reef.type),size = 3) + # add the point markers
    geom_polygon(data = hull.data, aes(x = MDS1, y = MDS2, fill = Reef.type, group = Reef.type), alpha = 0.40) + # add the convex hulls
    geom_text(data = species.scores, aes(x = NMDS1, y = NMDS2, label = category), 
              alpha=0.5, size = 4.5) +  
    labs(x = "NMDS1", y = "NDMS2") +
    scale_colour_manual(values = c('#225ea8','#a1dab4','#41b6c4','#edf8b1')) +
    scale_fill_manual(values = c('#225ea8','#a1dab4','#41b6c4','#edf8b1'))+
    coord_equal() +
    theme_diss() +
    theme(panel.border = element_rect(fill = NA),
          legend.title = element_blank(),
          legend.position = c(0.85,0.9)))

ggsave(nmds_plot_rt, filename = "outputs/nmds_rt.png", width = 10, height = 8)

## Reef type ANOSIM & SIMPER ----

group <- substr(data2$Reef.type, 1,56)
group[group=="1"] <- "Channel"
group[group=="2"] <- "Inner"
group[group=="3"] <- "Outer"
group[group=="4"] <- "Thila"

group.fac <- factor(group, levels = c("Channel", "Inner", "Outer", "Thila"))

# setting the seed means the 'random' outcome will be reproducible.
set.seed(123)
anosim(data1, data2$Reef.type, permutations = 999, distance = "bray", strata = NULL)

#ANOSIM statistic R: 0.2008 
#Significance: 0.001 

simper(data1, data2$Reef.type, permutations = 999, trace = FALSE)

## Atoll side NMDS ----

EWnmds_data <- all22sites_aggr %>% 
  filter(Reef.type %in% c("Channel", "Outer")) %>% 
  mutate(atoll.side = case_when(Site.ID %in% c("IM", "FK", "MBO", "MK", "FO") ~ "East",
                                Site.ID %in% c("VK", "MC", "MO", "FGO", "GO",
                                               "KO", "HC", "OC", "GC", "HW") ~ "West"))
# Reorder west before east
EWnmds_data$atoll.side <- factor(EWnmds_data$atoll.side, levels = c("West", "East"),
                                 labels = c("West", "East"))

# Subset dataframe, on which to base the ordination:
data1 <- EWnmds_data[,6:16]
# Descriptive columns
data2 <- EWnmds_data[,17]

set.seed(1506)
# Because the final result depends on the initial random placement of the points, 
# I`ll set a seed to make the results reproducible (get same outcome every time)

nmds2 <- metaMDS(data1, k = 2, trymax = 100, trace = F, autotransform = FALSE, distance="bray")
nmds2

stressplot(nmds2)
plot(nmds2)

ordiplot(nmds2, type = "n")
orditorp(nmds2, display = "species", col = "red", air = 0.01)
orditorp(nmds2, display = "site", cex = 1.1, air = 0.01)

data.scores2 <- as.data.frame(nmds2$points)
nmds_all2 <- cbind(data2, data.scores2)

species.scores <- as.data.frame(scores(nmds2, "species"))
species.scores$category <- rownames(species.scores)

# Make hulls
grp.a <- nmds_all2[nmds_all2$data2 == "West", ][chull(nmds_all2[nmds_all2$data2 == "West", c("MDS1", "MDS2")]), ]  
grp.b <- nmds_all2[nmds_all2$data2 == "East", ][chull(nmds_all2[nmds_all2$data2 == "East", c("MDS1", "MDS2")]), ]  
hull.data <- rbind(grp.a, grp.b)  #combine grp.a and grp.b
hull.data

# Plot with ggplot
(nmds_EW <- ggplot() + 
    geom_polygon(data = hull.data, aes(x = MDS1, y = MDS2, fill = data2, group = data2), alpha = 0.30) + # add the convex hulls
    geom_point(data = nmds_all2, aes(x = MDS1, y = MDS2, shape = data2, colour = data2),size = 3) + # add the point markers
    geom_text(data = species.scores, aes(x = NMDS1, y = NMDS2, label = category), 
              alpha=0.5, size = 4.5) +  # add the species labels
    #geom_text(data = nmds_all2, aes(x = NMDS1, y = NMDS2, label = survey.id), size = 1, vjust = 0) +  # add the site labels
    #annotate("text", x = -0.68, y = -0.76, label = "Ecklonia radiata", size = 3) + 
    labs(x = "NMDS1", y = "NDMS2") +
    scale_x_continuous(limits = c(-1.2, 0.8)) +
    scale_colour_manual(values = c("#a1dab4", "#2c7fb8")) +
    scale_fill_manual(values = c("#a1dab4", "#2c7fb8"))+
    coord_equal() +
    theme_diss() +
    theme(panel.border = element_rect(fill = NA),
          legend.title = element_blank(),
          #aspect.ratio = 1,
          legend.position = c(0.85,0.9)))

ggsave(nmds_EW, filename = "outputs/nmds_EW.png", width = 9, height = 8)

## Atoll side ANOSIM & SIMPER ----
group <- substr(nmds_all2$data2, 1,30)
group[group=="1"] <- "West"
group[group=="2"] <- "East"

group.fac <- factor(group, levels = c("West", "East"))

# setting the seed means the 'random' outcome will be reproducible.
set.seed(123)
anosim(data1, nmds_all2$data2, permutations = 999, distance = "bray", strata = NULL)

# ANOSIM statistic R: 0.2571 
# Significance: 0.005

simper(data1, nmds_all2$data2, permutations = 999, trace = FALSE)
# Bare rock

## Stacked barplots ----

### Year ----

# Need data in long format:
perc_site_long <- perc_site %>% 
  pivot_longer(cols = 6:16, names_to = "Cover", values_to = "Percent")

# Reorder categories
perc_site_long$Cover <- factor(perc_site_long$Cover,
                               levels = c("Other", "Bare.rock", "Rubble", "Sand",
                                          "Algae", "Sponge", "Soft.coral", "Dead.coral",
                                          "Diseased.coral", "Bleached.coral", "Hard.coral"),
                               labels = c("Other", "Bare rock", "Rubble", "Sand",
                                          "Algae", "Sponge", "Soft coral", "Dead coral",
                                          "Diseased coral", "Bleached coral", "Hard coral"))

(years_stack <- ggplot(perc_site_long, aes(x = year, y = Percent, fill = Cover)) +
    geom_bar(position = "fill", stat = "identity", width = 0.4) +
    scale_y_continuous(labels = scales::percent,
                       expand = c(0,0)) + 
    labs(y = "Benthic cover", x = "Year") +
    scale_colour_manual(values = c("Other" = "#9b9b7a","Bare rock" = "#baa587",
                                   "Rubble" = "#d9ae94", "Sand" = "#f1dca7",
                                   "Algae" = "#90be6d", "Sponge" = "#ffb700",
                                   "Soft coral" = "#d55d92", "Dead coral" = "#b1a7a6",
                                   "Diseased coral" = "#e2cfc4",
                                   "Bleached coral" = "#f4e1d9", "Hard coral" = '#f28c87')) +
    scale_fill_manual(values = c("Other" = "#9b9b7a","Bare rock" = "#baa587",
                                 "Rubble" = "#d9ae94", "Sand" = "#f1dca7",
                                 "Algae" = "#90be6d", "Sponge" = "#ffb700",
                                 "Soft coral" = "#d55d92", "Dead coral" = "#b1a7a6",
                                 "Diseased coral" = "#e2cfc4",
                                 "Bleached coral" = "#f4e1d9", "Hard coral" = '#f28c87')) +
    theme_diss() +
    theme(legend.position = "right",
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          #aspect.ratio = 1.5,
          panel.border = element_rect(fill = NA),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))
)

ggsave(years_stack, filename = "outputs/years_stack.png", width = 7, height = 6)

### Depth ----

# Use all 2022 data but same categories

# Need data in long format:
all22sites_aggr_long <- all22sites_aggr %>% 
  pivot_longer(cols = 6:16, names_to = "Cover", values_to = "Percent")

# Reorder depths
all22sites_aggr_long$Depth <- factor(all22sites_aggr_long$Depth,
                                     levels = c("5m", "10m", "15m"),
                                     labels = c("5m", "10m", "15m"))

# Reorder categories
all22sites_aggr_long$Cover <- factor(all22sites_aggr_long$Cover,
                                     levels = c("Other", "Bare.rock", "Rubble", "Sand",
                                                "Algae", "Sponge", "Soft.coral", "Dead.coral",
                                                "Diseased.coral", "Bleached.coral", "Hard.coral"),
                                     labels = c("Other", "Bare rock", "Rubble", "Sand",
                                                "Algae", "Sponge", "Soft coral", "Dead coral",
                                                "Diseased coral", "Bleached coral", "Hard coral"))

(depth_stack <- ggplot(all22sites_aggr_long, aes(x = Depth, y = Percent, fill = Cover)) +
    geom_bar(position = "fill", stat = "identity", width = 0.5) +
    scale_y_continuous(labels = scales::percent,
                       expand = c(0,0)) + 
    labs(y = "Benthic cover", x = "Depth") +
    scale_colour_manual(values = c("Other" = "#9b9b7a","Bare rock" = "#baa587",
                                   "Rubble" = "#d9ae94", "Sand" = "#f1dca7",
                                   "Algae" = "#90be6d", "Sponge" = "#ffb700",
                                   "Soft coral" = "#d55d92", "Dead coral" = "#b1a7a6",
                                   "Diseased coral" = "#e2cfc4",
                                   "Bleached coral" = "#f4e1d9", "Hard coral" = '#f28c87')) +
    scale_fill_manual(values = c("Other" = "#9b9b7a","Bare rock" = "#baa587",
                                 "Rubble" = "#d9ae94", "Sand" = "#f1dca7",
                                 "Algae" = "#90be6d", "Sponge" = "#ffb700",
                                 "Soft coral" = "#d55d92", "Dead coral" = "#b1a7a6",
                                 "Diseased coral" = "#e2cfc4",
                                 "Bleached coral" = "#f4e1d9", "Hard coral" = '#f28c87')) +
    theme_diss() +
    theme(legend.position = "right",
          legend.title = element_blank(),
          panel.border = element_rect(fill = NA),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
)

ggsave(depth_stack, filename = "outputs/depth_stack.png", width = 10, height = 8)

### Reef type ----

all22sites_aggr_long$Reef.type <- factor(all22sites_aggr_long$Reef.type,
                                         levels = c("Channel", "Outer", "Inner", "Thila"),
                                         labels = c("Channel", "Outer", "Inner", "Thila"))

(reeft_stack <- ggplot(all22sites_aggr_long, aes(x = Reef.type, y = Percent, fill = Cover)) +
    geom_bar(position = "fill", stat = "identity", width = 0.6) +
    scale_y_continuous(labels = scales::percent,
                       expand = c(0,0)) + 
    labs(y = "Benthic cover", x = "Reef type") +
    scale_colour_manual(values = c("Other" = "#9b9b7a","Bare rock" = "#baa587",
                                   "Rubble" = "#d9ae94", "Sand" = "#f1dca7",
                                   "Algae" = "#90be6d", "Sponge" = "#ffb700",
                                   "Soft coral" = "#d55d92", "Dead coral" = "#b1a7a6",
                                   "Diseased coral" = "#e2cfc4",
                                   "Bleached coral" = "#f4e1d9", "Hard coral" = '#f28c87')) +
    scale_fill_manual(values = c("Other" = "#9b9b7a","Bare rock" = "#baa587",
                                 "Rubble" = "#d9ae94", "Sand" = "#f1dca7",
                                 "Algae" = "#90be6d", "Sponge" = "#ffb700",
                                 "Soft coral" = "#d55d92", "Dead coral" = "#b1a7a6",
                                 "Diseased coral" = "#e2cfc4",
                                 "Bleached coral" = "#f4e1d9", "Hard coral" = '#f28c87')) +
    theme_diss() +
    theme(legend.position = "right",
          legend.title = element_blank(),
          panel.border = element_rect(fill = NA),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
)

ggsave(reeft_stack, filename = "outputs/reeft_stack.png", width = 10, height = 8)

### Atoll side ----

EWstack_data <- all22sites_aggr %>% 
  filter(Reef.type %in% c("Channel", "Outer")) %>% 
  mutate(atoll.side = case_when(Site.ID %in% c("IM", "FK", "MBO", "MK", "FO") ~ "East",
                                Site.ID %in% c("VK", "MC", "MO", "FGO", "GO",
                                               "KO", "HC", "OC", "GC", "HW") ~ "West"))
# Reorder west before east
EWstack_data$atoll.side <- factor(EWstack_data$atoll.side, levels = c("West", "East"),
                                  labels = c("West", "East"))

# Reshape
EWstack_data <- EWstack_data %>% 
  pivot_longer(cols = 6:16, names_to = "Cover", values_to = "Percent")

# Reorder categories
EWstack_data$Cover <- factor(EWstack_data$Cover,
                             levels = c("Other", "Bare.rock", "Rubble", "Sand",
                                        "Algae", "Sponge", "Soft.coral", "Dead.coral",
                                        "Diseased.coral", "Bleached.coral", "Hard.coral"),
                             labels = c("Other", "Bare rock", "Rubble", "Sand",
                                        "Algae", "Sponge", "Soft coral", "Dead coral",
                                        "Diseased coral", "Bleached coral", "Hard coral"))


(EW_stack <- ggplot(EWstack_data, aes(x = atoll.side, y = Percent, fill = Cover)) +
    geom_bar(position = "fill", stat = "identity", width = 0.4) +
    scale_y_continuous(labels = scales::percent,
                       expand = c(0,0)) + 
    labs(y = "Benthic cover", x = "Side of atoll") +
    scale_colour_manual(values = c("Other" = "#9b9b7a","Bare rock" = "#baa587",
                                   "Rubble" = "#d9ae94", "Sand" = "#f1dca7",
                                   "Algae" = "#90be6d", "Sponge" = "#ffb700",
                                   "Soft coral" = "#d55d92", "Dead coral" = "#b1a7a6",
                                   "Diseased coral" = "#e2cfc4",
                                   "Bleached coral" = "#f4e1d9", "Hard coral" = '#f28c87')) +
    scale_fill_manual(values = c("Other" = "#9b9b7a","Bare rock" = "#baa587",
                                 "Rubble" = "#d9ae94", "Sand" = "#f1dca7",
                                 "Algae" = "#90be6d", "Sponge" = "#ffb700",
                                 "Soft coral" = "#d55d92", "Dead coral" = "#b1a7a6",
                                 "Diseased coral" = "#e2cfc4",
                                 "Bleached coral" = "#f4e1d9", "Hard coral" = '#f28c87')) +
    theme_diss() +
    theme(legend.position = "right",
          legend.title = element_blank(),
          panel.border = element_rect(fill = NA),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
)

ggsave(EW_stack, filename = "outputs/EW_stack.png", width = 10, height = 8)

# Temporal analysis ----

temp_an_dep$Depth <- factor(temp_an_dep$Depth, levels = c("5m", "10m", "15m"),
                            labels = c("5m", "10m", "15m"))

temp_an_dep <- temp_an_dep %>% 
  mutate(across(c(Site, Site.ID, Reef.type, year), as.factor))

## Year LMM  ----

#with data aggr. by depth
lmm2 <- lmer(Hard.coral ~ year + (1|Site), data = temp_an_dep)
summary(lmm2)
plot(lmm2)  # no pattern!
simulateResiduals(lmm2, plot = T)  # fine!

# Calculate R2
r.squaredGLMM(lmm2)

# Plotting model predictions
pred.lmm2 <- ggpredict(lmm2, terms = c("year"))
plot(pred.lmm2)

## Year LMM + depth and reef type ----

lmm3 <- lmer(Hard.coral ~ year + Depth + Reef.type + year:Depth + year:Reef.type + (1|Site),
             data = temp_an_dep)
summary(lmm3)
plot(lmm3)
simulateResiduals(lmm3, plot = T) # no problems

# Get significance
anova(lmm3)

# Calculate R2
r.squaredGLMM(lmm3)

## Coral change over time plot ----

# Plotting change over time with raw data:
mean_cover <- temp_an %>% 
  group_by(year) %>% 
  summarise(Hard.coral.m = mean(Hard.coral), 
            Hard.coral.sd = sd(Hard.coral),
            n = n(),
            SE_c = sd(Hard.coral)/sqrt(n())) %>% 
  ungroup()

(time_change_cor <- ggplot(mean_cover, aes(x = year, y = Hard.coral.m)) +
    geom_point(size = 3.5, shape = 18) +
    scale_y_continuous(limits = c(0,45), expand = c(0,0)) +
    geom_errorbar(aes(ymin = Hard.coral.m - SE_c, 
                      ymax = Hard.coral.m + SE_c), colour = "black",
                  size = 0.5, width=0.05) +
    geom_point(data = perc_site, aes(x = year, y = Hard.coral, colour = year),
               position = position_jitter(width = 0.05), size = 2, alpha = 0.6) +
    scale_colour_manual(values = c("2019" = "#ffba49",
                                   "2021/22" = "#20a39e")) +
    scale_fill_manual(values = c("2019" = "#ffba49",
                                 "2021/22" = "#20a39e"))+
    labs(y = "Hard coral cover (%)") +
    theme_diss() +
    theme(axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.position = "none",
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
)

save_plot(time_change_cor, filename = "outputs/cor_change.png", base_height = 6, base_width = 6)

# Recovery dynamics ----

# Make column for positve/negative change
change <- change %>% 
  mutate(pos = perc.change >= 0)

## Perc. change bar plot ----

(cor.perc.change.plot <- ggplot(change, aes(x = reorder(Site, desc(perc.change)), 
                                            y = perc.change, fill = pos)) +
    geom_col(position = "identity", size = 0.25, alpha = 0.7) +
    geom_hline(yintercept = 0) +
    scale_fill_manual(values = c("#d42020","#1e8df4")) +
    #scale_colour_manual(values = c("#025ba4", "#d42020")) +
    scale_y_continuous(limits = c(-100, 450)) +
    labs(y = "% change in coral cover") +
    theme_diss() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
)

ggsave(cor.perc.change.plot, filename = "outputs/cor_perc_change.png", width = 9, height = 6)

## 2019 vs. % change ----

cor.ch.lm <- lm(perc.change ~ X2019.cover, data = change)
summary(cor.ch.lm)
plot(cor.ch.lm)
simulateResiduals(cor.ch.lm, plot = T)  # fine

cor.ch.lm.pred <- ggpredict(cor.ch.lm, terms = c("X2019.cover"))
plot(cor.ch.lm.pred)

(plot1 <- ggplot(cor.ch.lm.pred) + 
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                alpha = 0.2) +  # error band
    geom_point(data = change, size = 3, alpha = 0.5,                   
               aes(x = X2019.cover, y = perc.change, colour = pos)) + 
    scale_colour_manual(values = c("#d42020", "#1e8df4")) +
    annotate("text", x = 29, y = 390, label = "y = -8.88x + 236,", size = 5) +
    annotate("text", x = 29.5, y = 345,
             label = "paste(italic(R) ^ 2, \" = 0.36, p < 0.01\")", 
             parse = T, size = 5) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.6) +
    labs(x = "Coral cover in 2019 (%)", 
         y = "% change in coral cover") + 
    theme_diss() +
    theme(legend.position = "none",
          #axis.title.x = element_blank(),
          #plot.margin = unit(c(0.5, 0.1, 0.5, 0.8), "cm"),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
)

ggsave(plot1, filename = "outputs/perc.2019.png", width = 8, height = 6)

## 2019 vs. diff ----

cor.diff.lm <- lm(Difference ~ X2019.cover, data = change)
summary(cor.diff.lm)
plot(cor.diff.lm)
simulateResiduals(cor.diff.lm, plot = T)  # fine

cor.dif.lm.pred <- ggpredict(cor.diff.lm, terms = c("X2019.cover"))
plot(cor.dif.lm.pred)

(plot2 <- ggplot(cor.dif.lm.pred) + 
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                alpha = 0.2) +  # error band
    geom_point(data = change, size = 3, alpha = 0.5,                    
               aes(x = X2019.cover, y = Difference, colour = pos)) + 
    annotate("text", x = 29, y = 32, label = "y = -0.81x + 23.6,", size = 5) +
    annotate("text", x = 29.5, y = 27.5,
             label = "paste(italic(R) ^ 2, \" = 0.24, p < 0.05\")", 
             parse = T, size = 5) +
    scale_colour_manual(values = c("#1e8df4", "#d42020")) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.6) +
    labs(x = "Coral cover in 2019 (%)", 
         y = "Difference in coral cover") + 
    theme_diss() +
    theme(legend.position = "none",
          #plot.margin = unit(c(0.5, 0.1, 0.5, 0.8), "cm"),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
)

ggsave(plot2, filename = "outputs/diff.2019.png", width = 10, height = 8)

## 2022 vs. % change ----

cor.ch.lm2 <- lm(perc.change ~ X2022.cover, data = change)
summary(cor.ch.lm2)
plot(cor.ch.lm2)
simulateResiduals(cor.ch.lm2, plot = T)  # fine

cor.ch.lm2.pred <- ggpredict(cor.ch.lm2, terms = c("X2022.cover"))
plot(cor.ch.lm2.pred)

(plot3 <- ggplot(cor.ch.lm2.pred) + 
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                alpha = 0.2) +  # error band
    geom_point(data = change, size = 3,  alpha = 0.5,                   
               aes(x = X2022.cover, y = perc.change, colour = pos)) + 
    annotate("text", x = 20, y = 390, label = "y = 5.58x - 65.3,", size = 5) +
    annotate("text", x = 20.4, y = 345,
             label = "paste(italic(R) ^ 2, \" = 0.23, p < 0.05\")", 
             parse = T, size = 5) +
    scale_colour_manual(values = c("#1e8df4", "#d42020")) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.6) +
    labs(x = "Coral cover in 2021/22 (%)", 
         y = "% change in coral cover") + 
    theme_diss() +
    theme(legend.position = "none",
          #axis.title.y = element_blank(),
          #axis.title.x = element_blank(),
          # plot.margin = unit(c(0.5, 0.1, 0.5, 0.1), "cm"),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
)

ggsave(plot3, filename = "outputs/perc.2022.png", width = 10, height = 8)

## 2022 vs diff ----

cor.diff.lm2 <- lm(Difference ~ X2022.cover, data = change)
summary(cor.diff.lm2)
plot(cor.diff.lm2)
simulateResiduals(cor.diff.lm2, plot = T)  # fine

cor.dif.lm2.pred <- ggpredict(cor.diff.lm2, terms = c("X2022.cover"))
plot(cor.dif.lm2.pred)

(plot4 <- ggplot(cor.dif.lm2.pred) + 
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                alpha = 0.2) +  # error band
    geom_point(data = change, size = 3, alpha = 0.5,                 
               aes(x = X2022.cover, y = Difference, colour = pos)) + 
    annotate("text", x = 20, y = 32, label = "y = 0.89x - 14.0,", size = 5) +
    annotate("text", x = 20.8, y = 27.5,
             label = "paste(italic(R) ^ 2, \" = 0.55, p < 0.001\")", 
             parse = T, size = 5) +
    scale_colour_manual(values = c("#1e8df4", "#d42020")) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.6) +
    labs(x = "Coral cover in 2021/22 (%)", 
         y = "Difference in coral cover") + 
    theme_diss() +
    theme(legend.position = "none",
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
)

ggsave(plot4, filename = "outputs/diff.2022.png", width = 10, height = 8)

## Arrange all four plots in panel
(panel <- plot_grid(plot1, plot3, plot2, plot4, nrow = 2,
                    align = "vh", axis = "b",
                    labels = "AUTO", label_size = 18))
save_plot(panel, filename = "outputs/cor_change_final.png", base_height = 9, base_width = 12)

# Spatial analysis (2021/22) ----

all22sites$Depth <- factor(all22sites$Depth, levels = c("15m", "10m", "5m"),
                           labels = c("15m", "10m", "5m"))

## Depth on coral ----

dep1 <- lmer(Total.hard.coral ~ Depth + (1|Site.name), data = all22sites)
summary(dep1)
plot(dep1)

anova(dep1) # n.s.

## Reef type on coral ----

rt1 <- lmer(Total.hard.coral ~ Reef.type + (1|Site.name), data = all22sites)
summary(rt1)
plot(rt1)

anova(rt1)  # n.s.

## Atoll side on coral ----

# Filter out outer sites
EW_an <- all22sites %>% 
  filter(Reef.type %in% c("Channel", "Outer"))

# Make columnf or East or west
EW_an <- EW_an %>% 
  mutate(atoll.side = case_when(Site.ID %in% c("IM", "FK", "MBO", "MK", "FO") ~ "East",
                                Site.ID %in% c("VK", "MC", "MO", "FGO", "GO",
                                               "KO", "HC", "OC", "GC", "HW") ~ "West"))
# Reorder west before east
EW_an$atoll.side <- factor(EW_an$atoll.side, levels = c("West", "East"),
                           labels = c("West", "East"))

EW_an %>% group_by(atoll.side) %>%  summarise(length(unique(Site.ID)))

EW_an %>% group_by(atoll.side) %>% summarise(mean(Total.hard.coral))

ggplot(EW_an, aes(atoll.side, Total.hard.coral)) +
  geom_boxplot()
# Looks like eastern sites have higher coral cover

# Test statistically

EWmod <- lmer(Total.hard.coral ~ atoll.side + (1|Site.ID), data = EW_an)
summary(EWmod)  #n.s
plot(EWmod)

# Reef impacts ----

## Impacts over time ----

impacts.glm <- glm(Total.impacts ~ year, data = all_data, family = poisson)
summary(impacts.glm)
plot(impacts.glm)
simulateResiduals(impacts.glm, plot = T)  # v bad
testOverdispersion(impacts.glm)  # significant

1011.0/36 # 28.08333, so severe overdispersion

# Use quasipoisson
impacts.glm2 <- glm(Total.impacts ~ year, data = all_data, family = quasipoisson)
summary(impacts.glm2)
plot(impacts.glm2)

anova(impacts.glm2, test = "F")

## Impacts on coral cover ----

yr_im_lm <- lm(Hard.coral ~ Total.impacts + year, data = all_data)
summary(yr_im_lm)
plot(yr_im_lm)
simulateResiduals(yr_im_lm, plot = T)  # ok

## Impacts diff vs. coral change ----

# Is the increase/decrease in coral cover dependent on the change in impacts?

# Plot trend
plot(change$perc.change ~ change$impacts.diff) # negative trend

change.lm <- lm(perc.change ~ impacts.diff, data = change)
summary(change.lm)  # nearly significant
plot(change.lm)
simulateResiduals(change.lm, plot = T)  # ok (one outlier)

change.lm.pred <- ggpredict(change.lm, terms = c("impacts.diff"))
plot(change.lm.pred)

(coral_impacts_change <- ggplot(change.lm.pred) + 
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = change, size = 3, alpha = 0.5,                    
               aes(x = impacts.diff, y = perc.change, colour = pos)) + 
    annotate("text", x = 70, y = 380, label = "y = -1.18x + 117,", size = 6) +
    annotate("text", x = 70.9, y = 350,
             label = "paste(italic(R) ^ 2, \" = 0.15, p = 0.057\")", 
             parse = T, size = 6) +
    scale_colour_manual(values = c("#1e8df4", "#d42020")) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.6) +
    scale_x_continuous(breaks = c(-40, -20, 0, 20, 40, 60, 80, 100, 120, 140)) +
    labs(x = "Change in reef impacts since 2019", 
         y = "% change in coral cover") + 
    theme_diss() +
    theme(legend.position = "none",
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
)

ggsave(coral_impacts_change, filename = "outputs/cor_im_change_final.png", width = 10, height = 8)


## Impact type plots ----

# Combine 19 and 22 impacts data 

# Tidy 2019 impacts
impacts19 <- impacts19 %>% 
  select(- X, - X.1, - X.2) %>% 
  mutate(year = "2019")

# Tidy 2022 impacts & select only repeat sites
impacts22 <- impacts22 %>% 
  select(- Size, - Species, - Biomass, - Family, - Surveyors.Name, - Start.Time, - End.Time) %>% 
  filter(Site %in% c("OR", "OI", "FI", "LF", "FK", "HC", "MBO", "MBI", "GI",
                     "FO", "HW", "KO", "MC", "GO", "MI", "MO", "OC", "RDH", "HI")) %>% 
  rename(Replicate = Replicate.No.) %>% 
  mutate(year = "2021/22")

# Combine
all_impacts <- bind_rows(impacts19, impacts22)

# Count no of different impacts by site and year
all_impacts2 <- all_impacts %>% 
  group_by(Site, year, Reef.impacts) %>% 
  count(length(Reef.impacts))

# Put data into wide format to aggregate
allimpacts_wide <- all_impacts2 %>% 
  pivot_wider(names_from = Reef.impacts, values_from = n, values_fill = 0)

# Aggregate categories:
allimpacts_wide <- allimpacts_wide %>% 
  transmute(Site = Site, 
            'Drupella Damage' = `DRUPELLA DAMAGE`, 'Fishing lines' = `FISHING LINES`,
            'Trash' = TRASH, 'COTS/CS Damage' = `COTS/CS DAMAGE`, 'Anchor damage' = `ANCHOR DAMAGE`, 
            'Other damages' = `OTHER DAMAGES` + `SP TISSUE` + `PARTLY DEAD` + DEAD + PREDATION +
              `BLUE RINGS ON CORALS` + `WHITE BAND DISEASE` + `WHITE SPOT` 
            + `SPOTTED DISEASE` + DISEASE + `POSSIBLE DISEASE/ BACTERIA`
            + `BACTERIA DAMAGE` + `PINK DOTS` + `CORAL BLEACHING` + `SPOT BLEACHING`, 
            'Parrotfish damage' = `PARROT FISH DAMAGE` + `PARROTFISH DAMAGE`)

# Reshape to plot
allimpacts_long <- allimpacts_wide %>% 
  pivot_longer(cols = 3:9, names_to = "Impact type", values_to = "Impact count")

# Stacked barplot by year
(impacts_stack <- ggplot(allimpacts_long, aes(x = year, y = `Impact count`, 
                                              fill = `Impact type`)) +
    geom_bar(position = "fill", stat = "identity", width = 0.5) +
    scale_colour_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5')) +
    scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5')) +
    scale_y_continuous(labels = scales::percent,
                       expand = c(0,0)) + 
    labs(y = "Reef impacts") +
    theme_diss() +
    theme(legend.position = "right",
          axis.title.x = element_blank(),
          panel.border = element_rect(fill = NA))
)

save_plot(impacts_stack, filename = "outputs/impacts_stack.png", base_height = 8, base_width = 9)

# Sum up total no of impacts of different types per year
impact_count <- allimpacts_long %>% 
  group_by(year, `Impact type`) %>%  
  summarise(count = sum(`Impact count`))

(impacts_bar <- ggplot(impact_count, aes(x = reorder(`Impact type`,(desc(count))), y = count, 
                                         colour = year, fill = year )) +
    geom_bar(stat = "identity", position = position_dodge(0.6), width = 0.5, alpha = 0.8) +
    scale_colour_manual(values = c("#ffba49", "#20a39e")) +
    scale_fill_manual(values = c("#ffba49", "#20a39e")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0,205)) +
    labs(y = "Total no. of reef impacts", x ="\nType of impact") +
    theme_diss() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.85, 0.9),
          axis.text = element_text(size = 16),
          plot.margin = unit(c(0.5, 0.1, 0.5, 0.1), "cm")))

save_plot(impacts_bar, filename = "outputs/impacts_bar_new.png", base_height = 6)

# Plot total impacts seperately:

# Add collumn for total impacts
allimpacts_wide <- allimpacts_wide %>% 
  rowwise() %>% 
  mutate(`Total impacts` = sum(c_across(`Drupella Damage`:`Parrotfish damage`)))

# Reshape to plot
allimpacts_long <- allimpacts_wide %>% 
  pivot_longer(cols = 3:10, names_to = "Impact type", values_to = "Impact count")

# Sum up total no of impacts of different types per year
impact_count <- allimpacts_long %>% 
  group_by(year, `Impact type`) %>%  
  summarise(count = sum(`Impact count`))

total_impacts <- impact_count %>% 
  filter(`Impact type` == "Total impacts")

(total_impacts_bar <- ggplot(total_impacts, aes(x = `Impact type`, y = count, 
                                                colour = year, fill = year )) +
    geom_bar(stat = "identity", position = position_dodge(0.4), width = 0.3, alpha = 0.8) +
    scale_colour_manual(values = c("#ffba49", "#20a39e")) +
    scale_fill_manual(values = c("#ffba49", "#20a39e")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    labs(y = "Total no. of reef impacts\n", x ="\nType of impact") +
    theme_diss() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none",
          axis.text = element_text(size = 16),
          plot.margin = unit(c(0.5, 0.1, 0.5, 0.1), "cm")))

# Combine
(all_impacts <- plot_grid(impacts_bar, total_impacts_bar, ncol = 2, 
                          align = "hv", axis = "bt",
                          rel_widths = c(2, 1), 
                          labels = "AUTO", label_size = 18))
save_plot(all_impacts, filename = "outputs/all_impacts.png", base_height = 6, base_width = 10)

