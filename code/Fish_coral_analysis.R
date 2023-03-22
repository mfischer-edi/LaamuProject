# Laamu Coral Reef monitoring project
# Code for fish vs coral analysis

# Written by Mara Fischer
# E-mail: mf555@exeter.ac.uk
# 10/11/2022

# Packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(forcats)
library(lme4)
library(lmerTest)
library(DHARMa)
library(ggeffects)
library(MuMIn)

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

# Data wrangling ----

## Olhutholhu inside ----

# Fish data for this site was added late to the sharepoint,
# so this data wrangling was used to get fish abundance data
# from that site in order to combine it with the other data

# Get abundance data for Olhutholhu inside
OI <- read.csv("data/Olhutholhu_2022.csv", header = T)

# Create dataset with fish abundance
OI_abun <- OI %>% 
  group_by(Site, Depth, Replicate.No., Family, Species) %>% 
  summarise(count = length(Species))

write.csv(OI_abun, "data/OI_abun_raw.csv")

# Create one ignoring depth and rep
OI_abun_site <- OI %>% 
  group_by(Site, Family, Species) %>% 
  summarise(count = length(Species))

write.csv(OI_abun_site, "data/OI_abun_site.csv")

## 2022 data ----

fish22 <- read.csv("data/biomass_22.csv", header = T)

# Aggregate to get abundance and biomass in same data frame
fish_all22 <- fish22 %>% 
  group_by(site, depth, replicate_no, family, species, biomass) %>% 
  summarise(abundance = length(species)) %>% 
  ungroup()

## one row per species of a certain biomass,
## so if within the same replicate, the same species was observed
## at different sizes and hence has different biomasses, there are
## multiple rows for that certain species

sum(fish_all22$abundance)
sum(fish_all22$biomass)

## 2022 all raw fish data ----

raw_fish22 <- read.csv("data/raw_fish_2022.csv", header = T)

# Check data
unique(raw_fish22$Site) # all 29 sites
unique(raw_fish22$Reef.type) # the 4 reef types
unique(raw_fish22$Depth)  # the 3 depths
unique(raw_fish22$Replicate.No.) # 2 reps
unique(raw_fish22$Family)  # 9 families
unique(raw_fish22$Species) # 106 species

# Aggregate to get abundance and biomass in same data frame
agg_fish22 <- raw_fish22 %>% 
  group_by(Site, Depth, Replicate.No., Family, Species, Biomass) %>% 
  summarise(abundance = length(Species),
            Site.name = unique(Site.name)) %>% 
  ungroup()
## now have the abundance of a species of a given biomass at each rep

# Create dataset with fish abundance without biomass
abun_fish22 <- raw_fish22 %>% 
  group_by(Site, Depth, Replicate.No., Family, Species) %>% 
  summarise(abundance = length(Species),
            Site.name = unique(Site.name),
            reef.type = unique(Reef.type)) %>% 
  ungroup()

write.csv(abun_fish22, "data/abun_rep_sp_2022.csv")

# Create one ignoring depth and rep
abun_site22 <- raw_fish22 %>% 
  group_by(Site, Family, Species) %>% 
  summarise(count = length(Species),
            Site.name = unique(Site.name),
            reef.type = unique(Reef.type),
            year = "2021/22") %>% 
  ungroup()

write.csv(abun_site22, "data/abun_site_sp_2022.csv")

# Dataset for overall abundance & biomass (ignoring species)
abun_site_agg22 <- raw_fish22 %>% 
  group_by(Site, Depth, Replicate.No.) %>% 
  summarise(abundance = length(Species),
            Site.name = unique(Site.name),
            reef.type = unique(Reef.type),
            year = "2021/22") %>% 
  ungroup()

write.csv(abun_site_agg22, "data/abun_site_2022.csv")

## Biomass data ----

# Select the six families we have biomass data for

biom_data <- raw_fish22 %>% 
  filter(Family %in% c("SNAP", "EMP", "GRO", "RAB", "T&J", "PAR"))

# Remove NAs
biom_data <- na.omit(biom_data)

# Add status column
biom_data <- biom_data %>% 
  mutate(Status = case_when(Family == "PAR" ~ "ECO", Family == "EMP" | 
                            Family == "GRO" | Family == "SNAP" |
                            Family == "T&J" | Family == "RAB" ~ "COM"),
         Rep = case_when(Replicate.No. == 1 ~ "Rep1",
                         Replicate.No. == 2 ~ "Rep2"))

# Aggreagte at site level (for GIS)
biom_site <- biom_data %>% 
  group_by(Site.name) %>% 
  summarise(Reef.type = unique(Reef.type),
            Biomass = mean(Biomass)) %>% 
  ungroup()

#write.csv(biom_site, "data/biomass_by_site.csv")

# Aggregate sum rather than average
biom_site_sum <- biom_data %>% 
  group_by(Site.name) %>% 
  summarise(Reef.type = unique(Reef.type),
            Biomass = sum(Biomass)/1000) %>% 
  ungroup()

## gives the total site biomass (per 2000m2) in kgs
## most match up with Alicia's table, but some don't (especially higher ones)

# Check for Fushi kandu
fushi <- raw_fish22 %>% 
  filter(Site == "FK")

fushi <- na.omit(fushi)

fushi2 <- fushi %>% 
  group_by(Site.name) %>% 
  summarise(Biomass = sum(Biomass)/1000) %>% 
  ungroup()
## 284.99 which matches with Alicia's so she left in all families, including SL
## whereas I only included the six

unique(fushi$Family) # yep included SL
  
# Check for Laamafaruhaa
lfh <- raw_fish22 %>% 
  filter(Site == "LF")

unique(lfh$Family) # inlcudes Trig and But (but they don't have biomass??)

lfh <- na.omit(lfh)

lfh %>%  summarise(sum(Biomass)/1000)

# Aggregate at rep level to combine with coral
biom_rep <- biom_data %>% 
  group_by(Site.name, Depth, Rep) %>% 
  summarise(Reef.type = unique(Reef.type),
            Biomass = mean(Biomass)) %>% 
  ungroup()

# now have mean biomass at replicate level (g per 500m2)

# Calculate the total biomass at each rep rather than the average
biom_rep_sum <- biom_data %>% 
  group_by(Site.name, Depth, Rep) %>% 
  summarise(Reef.type = unique(Reef.type),
            Biomass = sum(Biomass)) %>% 
  ungroup()

# Divide into eco and com

eco_biom <- biom_data %>% 
  filter(Family == "PAR")

eco_rep <- eco_biom %>% 
  group_by(Site.name, Depth, Rep) %>% 
  summarise(Reef.type = unique(Reef.type),
            Biomass = mean(Biomass)) %>% 
  ungroup()

com_biom <- biom_data %>% 
  filter(Family %in% c("SNAP", "EMP", "GRO", "RAB", "T&J"))

com_rep <- com_biom %>% 
  group_by(Site.name, Depth, Rep) %>% 
  summarise(Reef.type = unique(Reef.type),
            Biomass = mean(Biomass)) %>% 
  ungroup()

## Import 22 coral data

coral_all22 <- read.csv("data/perc_all_final.csv", header = T) # full 2022 raw data (28 sites)

# Factor depth
coral_all22$Depth <- factor(coral_all22$Depth, levels = c("5m", "10m", "15m"),
                           labels = c("5m", "10m", "15m"))

# Add column for total cover & Remove unneeded columns & get column for rep
coral_all22 <- coral_all22 %>% 
  mutate(Total.coral = Acropora + Coral) %>% 
  select(- Tape, - Wand, - Shadow, - Sum.excl.TWS) %>% 
  separate(Photo.Name, into = c("S", "D", "Rep", "Photo"), sep = "_")

# Check Rep column:
unique(coral_all22$Rep) # 2 reps

coral_all22$Rep <- as.factor(coral_all22$Rep)

# Remove unneeded columns
coral_all22 <- coral_all22 %>% 
  select(- S, - D, - Photo)

## Aggregate mean coral at site level
site_coral <- coral_all22 %>% 
  group_by(Site.name) %>% 
  summarise(Reef.type = unique(Reef.type), 
            Total.coral = mean(Total.coral)) %>%  
  ungroup()

#write.csv(site_coral, "coral22_site.csv")

## Aggregate mean coral at replicate level
coral_rep <- coral_all22 %>% 
  group_by(Site.name, Depth, Rep) %>% 
  summarise(Reef.type = unique(Reef.type), 
            Total.coral = mean(Total.coral)) %>%  
  ungroup()

# Add new depth column without m in it
coral_rep <- coral_rep %>% 
  mutate(Depth2 = case_when(Depth == "5m" ~ 5,
                            Depth == "10m" ~ 10,
                            Depth == "15m" ~ 15))

coral_rep <- coral_rep %>% 
  select(- Depth) %>% 
  rename(Depth = Depth2)

unique(coral_rep$Site.name)
unique(biom_rep$Site.name)

# Remove pink thila as no coral data for that
biom_rep <- biom_rep %>% 
  filter(Site.name != "Pink thila")

# Combine biomass and coral at rep level data *for 28 sites:

biom_coral <- left_join(biom_rep, coral_rep, by = c("Site.name","Depth","Rep"))

biom_coral <- biom_coral %>% 
  select(-Reef.type.x) %>% 
  rename(Reef.type = Reef.type.y)

write.csv(biom_coral, "data/biom_coral.csv")

# Combine total biomass and coral at rep level
biom_sum_coral <- left_join(biom_rep_sum, coral_rep, by = c("Site.name","Depth","Rep"))

# convert biomass to kg
biom_sum_coral <- biom_sum_coral %>% 
  mutate(Biomass.kg = Biomass/1000)

# combine eco_rep and coral
eco_biom_coral <- left_join(eco_rep, coral_rep, by = c("Site.name","Depth","Rep"))

eco_biom_coral <- eco_biom_coral %>% 
  filter(Site.name != "Pink thila") %>% 
  select(-Reef.type.x) %>% 
  rename(Reef.type = Reef.type.y)

write.csv(eco_biom_coral, "data/eco_biom_coral.csv")

# combine com rep and coral
com_biom_coral <- left_join(com_rep, coral_rep, by = c("Site.name","Depth","Rep"))

com_biom_coral <- com_biom_coral %>% 
  filter(Site.name != "Pink thila") %>% 
  select(-Reef.type.x) %>% 
  rename(Reef.type = Reef.type.y)

write.csv(com_biom_coral, "data/com_biom_coral.csv")

## 2019 abundance data ----

abun19 <- read.csv("data/abun_site_species_2019.csv", header = T)

# full abundance data for 2019 but aggregated by site and species,
# no distinction between depths or replicates

unique(abun19$site_id) # 20 repeat sites
unique(abun19$reef_type) # 4
unique(abun19$family) # 7, but including SL
unique(abun19$species) # 96 spp.

## Combine 2019 & 2022 abundance data ----

# Rename columns so they are the same
abun_site22 <- abun_site22 %>% 
  rename(Site.id = Site)

abun19 <- abun19 %>% 
  rename(Site.id = site_id, Site.name = site, reef.type = reef_type,
         Family = family, Species = species)

# Select only the 20 repeat sites from 2022 data
abun_site22_repeat <- abun_site22 %>% 
  filter(Site.id %in% c("OR", "OI", "FI", "LF", "FK", "HC", "MBO", "MBI", "GI", "PT",
                      "FO", "HW", "KO", "MC", "GO", "MI", "MO", "OC", "RDH", "HI")) 
  
# Change into factors
abun_site22_repeat <- abun_site22_repeat %>% 
  mutate(across(c(Site.id, Site.name, reef.type, year), as.factor))

abun19 <- abun19 %>% 
  mutate(across(c(Site.id, Site.name, reef.type, year), as.factor))

# Combine both years
abun_19_22 <- bind_rows(abun19, abun_site22_repeat)

unique(abun_19_22$Family)

# Select the six families that were surveyed in both years (and remove SL)
abun_19_22 <- abun_19_22 %>% 
  filter(Family %in% c("BUT", "EMP", "GRO", "PAR", "SNAP", "TRIG"))

# Add column for commercially vs ecologically important
abun_19_22 <- abun_19_22 %>% 
  mutate(Status = case_when(Family == "BUT" | Family == "PAR" | Family == "TRIG" ~ "ECO",
                            Family == "EMP" | Family == "GRO" | Family == "SNAP" ~ "COM"))

# Save this as final abundance dataset!

write.csv(abun_19_22, "data/FISH_ABUN_19_22.csv")

#abun_19_22 <- read.csv("data/FISH_ABUN_19_22.csv", header = T)

# Make two separate datasets for ECO and COM fish

eco_abun <- abun_19_22 %>% 
  filter(Status == "ECO")

write.csv(eco_abun, "data/eco_abun.csv")

com_abun <- abun_19_22 %>% 
  filter(Status == "COM")

write.csv(com_abun, "data/com_abun.csv")

# Combine coral and fish data ----

# Import coral data aggregated by site

coral <- read.csv("data/perc_site.csv", header = T)

coral <- coral %>% 
  rename(Site.id = Site.ID, Site.name = Site, reef.type = Reef.type)

# Join to fish abun

#abun_coral <- dplyr::left_join(abun_19_22, coral, by = "Site.id")
## this doesn't really work

# Aggregate total abundance at site level (all species combined)

abun_19_22_site <- abun_19_22 %>% 
  group_by(Site.id, year) %>% 
  summarise(Site.name = unique(Site.name), reef.type = unique(reef.type),
            abundance = sum(count)) %>% 
  ungroup()

# Double check this:
#test <- abun_19_22 %>% 
 # filter(Site.id == "FI")

#sum(test$count)
## adds up yes

# Remove pink thila as no coral data for that
abun_19_22_site <- abun_19_22_site %>% 
  filter(Site.id != "PT")

# Simplify coral dataset
coral2 <- coral %>% 
  select(Site.id, year, Hard.coral)

# Add coral data to this instead

abun_coral_site <- left_join(abun_19_22_site, coral2, by = c("Site.id", "year"))

## Dataset with two rows per site, one for each year

write.csv(abun_coral_site, "data/abun_coral_site.csv")

## Join coral to eco and com datasets

# Aggregate total abundance at site level (for eco and com)

eco_abun_site <- eco_abun %>% 
  group_by(Site.id, year) %>% 
  summarise(Site.name = unique(Site.name), reef.type = unique(reef.type),
            abundance = sum(count)) %>% 
  ungroup()

com_abun_site <- com_abun %>% 
  group_by(Site.id, year) %>% 
  summarise(Site.name = unique(Site.name), reef.type = unique(reef.type),
            abundance = sum(count)) %>% 
  ungroup()

# Remove pink thila as no coral data for that
eco_abun_site <- eco_abun_site %>% 
  filter(Site.id != "PT")

com_abun_site <- com_abun_site %>% 
  filter(Site.id != "PT")

# Combine with coral

eco_abun_coral <- left_join(eco_abun_site, coral2, by = c("Site.id", "year"))
write.csv(eco_abun_coral, "data/eco_abun_coral.csv")
com_abun_coral <- left_join(com_abun_site, coral2, by = c("Site.id", "year"))
write.csv(com_abun_coral, "data/com_abun_coral.csv")

# ANALYSIS ----

## Coral vs abun (all) ----

abun_coral_site <- read.csv("data/abun_coral_site.csv", header = T)

plot(abun_coral_site$Hard.coral ~ abun_coral_site$abundance)

# one big outlier for abundance (maavah outside)

# Remove outlier and plot again
abun_coral_site2 <- abun_coral_site %>% 
  filter(Site.id != "MO")

plot(abun_coral_site2$Hard.coral ~ abun_coral_site2$abundance) # weak pos trend

# Simple LM

lm1 <- lm(Hard.coral ~ abundance, data = abun_coral_site2)
summary(lm1)  # abundance has a significant weak positive effect on coral cover
plot(lm1)
# Adj. R2: 0.1008

## abundance alone only explains 10% of the variation in coral cover

# Incorporate temporal aspect
lm2 <- lm(Hard.coral ~ abundance + year + abundance*year, data = abun_coral_site2)
summary(lm2)

drop1(lm2, test = "F")

# non-sig. interaction terms indicates that the effect of abundance 
# on coral cover does not differ between the years

# Include year without interaction
lm3 <- lm(Hard.coral ~ abundance + year, data = abun_coral_site2)
summary(lm3)
plot(lm3)
# Adj. R2: 0.359

## Including year as a fixed effect helps explain more variation in coral cover
## since coral cover differs significantly over time, so abundance and year
## together explain 36% of the variation in coral cover

# Plot lm1 (without year)
(abun.p <- ggplot(abun_coral_site2, aes(x = abundance, y = Hard.coral)) +
    geom_point(size = 3, alpha = 0.6) +                                
    labs(x = expression(paste('Fish abundance (individuals per 2000' , m^2,')')),
         y = "Coral cover (%)",
         title = "A") +
    stat_smooth(method = "lm", colour = "black") +  
    scale_y_continuous(limits = c(0,45)) +
    theme_diss() +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
)

ggsave2(abun.p, filename = "outputs/abun_coral_lm.png", width = 7, height = 5)

# Plot lm2 (with random slopes for year)
(abun.year <- ggplot(abun_coral_site2, aes(x = abundance, y = Hard.coral)) +
    geom_point(aes(colour = year)) +                                
    labs(x = "Fish abundance (ind. per 2000 m2)", y = "Coral cover (%)") +
    stat_smooth(method = "lm", aes(fill = year, colour = year)) +   
    scale_colour_manual(values = c("#FFC125", "#36648B")) +
    scale_fill_manual(values = c("#FFC125", "#36648B")) +
    theme_diss() +
    theme(legend.title = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))))

ggsave2(abun.year, filename = "outputs/abun__year_coral_lm.png", width = 10, height = 7)

## Coral vs abun (ECO/COM) ----

# Load data divided into com and eco
eco_abun_coral <- read.csv("data/eco_abun_coral.csv", header = T)
com_abun_coral <- read.csv("data/com_abun_coral.csv", header = T)

# Remove Maavah outside outlier
eco_abun_coral2 <- eco_abun_coral %>% 
  filter(Site.id != "MO")

com_abun_coral2 <- com_abun_coral %>% 
  filter(Site.id != "MO")

### ECO

eco1 <- lm(Hard.coral ~ abundance, data = eco_abun_coral2)
summary(eco1) # weak sig., v low R2
plot(eco1)

eco1_pred <- ggpredict(eco1, terms = c("abundance[all]"))
plot(eco1_pred)

eco2 <- lm(Hard.coral ~ abundance + year, data = eco_abun_coral2)
summary(eco2) # abun no longer sig., higher R2

eco3 <- lm(Hard.coral ~ abundance*year, data = eco_abun_coral2)
summary(eco3) # n.s. interaction

(eco.abun <- ggplot(eco_abun_coral, aes(x = abundance, y = Hard.coral)) +
    geom_point(size = 2) +                                
    labs(x = "Abundance of ecologically important species (ind. per 2000 m2)",
         y = "Coral cover (%)") +
    stat_smooth(method = "lm", colour = "black") +    
    theme_diss() +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
)

ggsave2(eco.abun, filename = "outputs/eco_abun_lm.png", width = 10, height = 7)

### COM

com1 <- lm(Hard.coral ~ abundance, data = com_abun_coral2)
summary(com1) # not sig.
plot(com1)

com1_pred <- ggpredict(com1, terms = c("abundance[all]"))
plot(com1_pred)

com2 <- lm(Hard.coral ~ abundance + year, data = com_abun_coral2)
summary(com2) # year sig., much higher R2

com3 <- lm(Hard.coral ~ abundance*year, data = com_abun_coral2)
summary(com3) # n.s. interaction

(com.abun <- ggplot(com_abun_coral2, aes(x = abundance, y = Hard.coral)) +
    geom_point(size = 2) +                                
    labs(x = "Abundance of commercially important species (ind. per 2000 m2)",
         y = "Coral cover (%)") +
    stat_smooth(method = "lm", colour = "black") +    
    theme_diss() +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
)

ggsave2(com.abun, filename = "outputs/com_abun_lm.png", width = 10, height = 7)

# Arrange in panel
library(gridExtra)
panel <- grid.arrange(eco.abun, com.abun, nrow = 2)
ggsave(panel, filename = "outputs/eco_com_abun.png", width = 10, height = 14)

# Plot eco and com in same graph
cols <- c("eco"="#b80c09","com"="#0b4f6c")

(eco.com.abun <- ggplot() +
    geom_line(data = eco1_pred, aes(x = x, y = predicted, colour = "eco"), ) +         
    geom_ribbon(data = eco1_pred, alpha = 0.3,
                aes(x = x, ymin = predicted - std.error, 
                    ymax = predicted + std.error, fill = "eco")) + 
    geom_point(data = eco_abun_coral, aes(x = abundance, y = Hard.coral, colour = "eco"), 
               size = 3, alpha = 0.6) + 
    geom_line(data = com1_pred, aes(x = x, y = predicted, colour = "com")) +         
    geom_ribbon(data = com1_pred, alpha = 0.3,
                aes(x = x, ymin = predicted - std.error, 
                                      ymax = predicted + std.error,
                                      fill = "com")) + 
    geom_point(data = com_abun_coral2, aes(x = abundance, y = Hard.coral, colour = "com"), 
               size = 3, alpha = 0.6) + 
    scale_y_continuous(limits = c(0,45)) +
    scale_colour_manual(name = "Status", values = cols, guide = guide_legend(fill = NULL,colour = NULL)) + 
    scale_fill_manual(name = "Stauts", values = cols, guide="none") +
    labs(x = expression(paste('Fish abundance (individuals per 2000' , m^2,')')),
         y = "Coral cover (%)",
         title = "B") +
    theme_diss() +
    theme(legend.position = c(0.9,1),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
)

ggsave2(eco.com.abun, filename = "outputs/eco_com_abun.png", height = 5, width = 7)

## Coral vs biomass (all) ----

biom_coral <- read.csv("data/biom_coral.csv", header = T)

plot(biom_coral$Total.coral ~ biom_coral$Biomass)

bio1 <- lmer(Total.coral ~ Biomass + (1|Site.name), data = biom_coral)
summary(bio1) # sig.
plot(bio1) # looks okay
bio1_sim <- simulateResiduals(bio1)  # problem with residuals vs predicted
plot(bio1_sim)
bio1_sim2 <- recalculateResiduals(bio1_sim, group = biom_coral$Site.name)
plot(bio1_sim2)  # fine now

r.squaredGLMM(bio1)
# R2m        R2c
# 0.04380573 0.4033918

bio1_pred <- ggpredict(bio1, terms = c("Biomass[all]"))
plot(bio1_pred)

(biom_cor_plot <- ggplot() +
    geom_line(data = bio1_pred, aes(x = x, y = predicted)) +         
    geom_ribbon(data = bio1_pred, aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  
    geom_point(data = biom_coral, aes(x = Biomass, y = Total.coral),
               size = 3, alpha = 0.6) +
    #scale_x_continuous(expand = c(0.1, 0.1)) + 
    #scale_y_continuous(expand = c(0.1, 0.1)) +
    labs(title = "A",
         y = "Coral cover (%)",
         x = expression(paste('Fish biomass (g/500 ' , m^2,')'))) +
    theme_diss() +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
)

ggsave2(biom_cor_plot, filename = "outputs/biom_cor_plot.png", width = 7, height = 5)

## Coral vs biomass (using sums) ----

plot(biom_sum_coral$Total.coral ~ biom_sum_coral$Biomass.kg)

bio.sum1 <- lmer(Total.coral ~ Biomass.kg + (1|Site.name), data = biom_sum_coral)
summary(bio.sum1) # not sig.
plot(bio.sum1) # looks okay
bio.sum1.sim <- simulateResiduals(bio.sum1, plot = T)  # good

r.squaredGLMM(bio.sum1)
#    R2m      R2c
# 0.001607897 0.408947

bio.sum1_pred <- ggpredict(bio.sum1, terms = c("Biomass.kg[all]"))
plot(bio.sum1_pred)

(biom_sum_cor_plot <- ggplot() +
    geom_line(data = bio.sum1_pred, aes(x = x, y = predicted)) +         
    geom_ribbon(data = bio.sum1_pred, aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  
    geom_point(data = biom_sum_coral, aes(x = Biomass.kg, y = Total.coral),
               size = 3, alpha = 0.6) +
    #scale_x_continuous(expand = c(0.1, 0.1)) + 
    #scale_y_continuous(expand = c(0.1, 0.1)) +
    labs(title = "A",
         y = "Coral cover (%)",
         x = expression(paste('Fish biomass (kg 500 ' , m^2,')'))) +
    theme_diss() +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
)

ggsave2(biom_sum_cor_plot, filename = "outputs/biom_sum_cor_plot.png", width = 7, height = 5)

## Coral vs biomass (ECO/COM) ----

### ECO 

eco_biom_coral <- read.csv("data/eco_biom_coral.csv", header = T)

eco_bio_m <- lmer(Total.coral ~ Biomass + (1|Site.name), data = eco_biom_coral)
summary(eco_bio_m) # just sig.
plot(eco_bio_m) # okay
eco_bio_sim <- simulateResiduals(eco_bio_m)  # fine

anova(eco_bio_m)

r.squaredGLMM(eco_bio_m)
# R2m        R2c
# 0.05044731 0.3838314

eco_bio_m_pred <- ggpredict(eco_bio_m, terms = c("Biomass[all]"))
plot(eco_bio_m_pred)

### COM 

com_biom_coral <- read.csv("data/com_biom_coral.csv", header = T)

com_bio_m <- lmer(Total.coral ~ Biomass + (1|Site.name), data = com_biom_coral)
summary(com_bio_m) # not sig. 
plot(com_bio_m) # okay
com_bio_sim <- simulateResiduals(com_bio_m)  # fine
plot(com_bio_sim) # problem with residuals vs pred
com_bio_sim <- recalculateResiduals(com_bio_sim, group = com_biom_coral$Site.name)
plot(com_bio_sim) # still problem

r.squaredGLMM(com_bio_m)

com_bio_m_pred <- ggpredict(com_bio_m, terms = c("Biomass[all]"))
plot(com_bio_m_pred)

# Plot eco and com in same graph
cols <- c("eco"="#b80c09","com"="#0b4f6c")

# Plot both together in same plot
(eco_com_plot <- ggplot() +
    geom_line(data = eco_bio_m_pred, aes(x = x, y = predicted, colour = "eco")) +         
    geom_ribbon(data = eco_bio_m_pred, alpha = 0.3,
                aes(fill = "eco", x = x, ymin = predicted - std.error, 
                    ymax = predicted + std.error)) + 
    geom_point(data = eco_biom_coral, aes(x = Biomass, y = Total.coral, colour = "eco"),
               size = 3, alpha = 0.6) +
    geom_line(data = com_bio_m_pred, aes(x = x, y = predicted, colour = "com")) +         
    geom_ribbon(data = com_bio_m_pred, alpha = 0.3,
                aes(x = x, ymin = predicted - std.error, 
                    ymax = predicted + std.error, fill = "com")) + 
    geom_point(data = com_biom_coral, aes(x = Biomass, y = Total.coral, colour = "com"),
               size = 3, alpha = 0.6) +
    scale_colour_manual(name = "Status", values = cols, guide = guide_legend(fill = NULL,colour = NULL)) + 
    scale_fill_manual(name = "Status", values = cols, guide="none") +
    labs(title = "B",
         y = "Coral cover (%)",
         x = expression(paste('Fish biomass (g/500 ' , m^2,')'))) +
    theme_diss() +
    theme(legend.position = c(.9,.9),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
)

ggsave2(eco_com_plot, filename = "outputs/eco_com_biom.png", width = 7, height = 5)

