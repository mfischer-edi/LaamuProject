# Laamu Coral Reef monitoring project
# Code for coral community composition analysis

# Written by Mara Fischer
# E-mail: mf555@exeter.ac.uk
# 20/10/2022

# Packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(forcats)
library(lme4)
library(lmerTest)
library(DHARMa)

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

# Helper functions ----

## SummarySE function
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
## data: a data frame.
## measurevar: the name of a column that contains the variable to be summariezed
## groupvars: a vector containing names of columns that contain grouping variables
## na.rm: a boolean that indicates whether to ignore NA's
## conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Load data ----

perc_all22 <- read.csv("data/perc_all_final.csv", header = T) # full 2022 raw data (28 sites)
all22sites <- read.csv("data/all22sites.csv", header = T)  # full data, aggr. by depth

# Data exploration with raw data ----

length(unique(perc_all22$Site.name))  # 28 sites

# Mean Acropora cover: 
mean(perc_all22$Acropora)  # 4.889018

# Mean other coral cover:
mean(perc_all22$Coral)  # 22.60464

# Mean total hard coral cover:
mean(perc_all22$Acropora) + mean(perc_all22$Coral) # 27.49366

# Percentage of acropora out of total:
4.889018/27.49366 * 100  # 17.78235 %

# Percentage of non-acropora out of total:
22.60464/27.49366 * 100  # 82.21765 %

# Data manipulation raw data ----

# Factor depth
perc_all22$Depth <- factor(perc_all22$Depth, levels = c("5m", "10m", "15m"),
                           labels = c("5m", "10m", "15m"))

# Add column for total cover & Remove unneeded columns & get column for rep
perc_all22 <- perc_all22 %>% 
  mutate(Total.coral = Acropora + Coral) %>% 
  select(- Tape, - Wand, - Shadow, - Sum.excl.TWS) %>% 
  separate(Photo.Name, into = c("S", "D", "Rep", "Photo"), sep = "_")

# Check Rep column:
unique(perc_all22$Rep) # 2 reps

# Remove unneeded columns
perc_all22 <- perc_all22 %>% 
  select(- S, - D, - Photo)

## Aggregate at site level (for GIS) ----

perc_site <- perc_all22 %>% 
  group_by(Site.name) %>% 
  summarise(Reef.type = unique(Reef.type), 
            Acropora = mean(Acropora), Other.coral = mean(Coral),
            Total.coral = mean(Total.coral), 
            Bleached.Acropora = mean(Bleached.Acropora),
            Bleached.other.coral = mean(Bleached.Coral),
            Total.bleached = mean(Bleached.Acropora) + mean(Bleached.Coral),
            Acropora.perc = Acropora/Total.coral * 100) %>% 
  ungroup()

write.csv(perc_site, "data/Acro_site.csv")

## Aggregate at replicate level ----
perc_rep <- perc_all22 %>% 
  group_by(Site.name, Depth, Rep) %>% 
  summarise(Reef.type = unique(Reef.type), 
            Acropora = mean(Acropora), Other.coral = mean(Coral),
            Total.coral = mean(Total.coral), 
            Bleached.Acropora = mean(Bleached.Acropora),
            Bleached.other.coral = mean(Bleached.Coral),
            Total.bleached = mean(Bleached.Acropora) + mean(Bleached.Coral)) %>% 
  ungroup()

# Add column for total coral incl. bleached

perc_rep <- perc_rep %>% 
  mutate(Total.incl.bleached = Total.coral + Total.bleached)

# Check averages still match up
mean(perc_rep$Acropora)
mean(perc_rep$Other.coral)
mean(perc_rep$Total.coral)
## all the same as before!

# Laamu-wide averages for bleached coral
mean(perc_rep$Bleached.Acropora)  # 0.695625
mean(perc_rep$Bleached.other.coral)  # 0.07392857
mean(perc_rep$Total.bleached)  # 0.7695536

# Percentage of acropora out of total bleached:
0.695625/0.7695536 * 100  # 90.39331 %

# Percentage of non-acropora out of total bleached:
0.07392857/0.7695536 * 100   # 9.606682 %

# Laamu-wide average for total coral incl. bleached
mean(perc_rep$Total.incl.bleached)  # 28.26321

# Percentage of bleached coral out of total:
0.7695536/28.26321 * 100   # 2.72281 %

# Percentage of live coral out of total:
27.49366/28.26321 * 100   # 97.2772 %

## Calculate mean cover + SE ----

# From raw data:
# Calculate mean, sd and SE for both Acropora and other coral
mean_cover <- perc_all22 %>% 
  summarise(Acropora.mean = mean(Acropora), 
            Acropora.sd = sd(Acropora),
            n = n(),
            Acropora.SE = sd(Acropora)/sqrt(n()),
            Coral.mean = mean(Coral), 
            Coral.sd = sd(Coral),
            Coral.SE = sd(Coral)/sqrt(n()),
            Total.mean = mean(Total.coral),
            Total.sd = sd(Total.coral),
            Total.SE = sd(Total.coral)/sqrt(n)) %>% 
  ungroup()

# From data aggregated by replicate:
mean_cover2 <- perc_rep %>% 
  summarise(Acropora.mean = mean(Acropora), 
            Acropora.sd = sd(Acropora),
            n = n(),
            Acropora.SE = sd(Acropora)/sqrt(n()),
            Coral.mean = mean(Other.coral), 
            Coral.sd = sd(Other.coral),
            Coral.SE = sd(Other.coral)/sqrt(n()),
            Total.mean = mean(Total.coral),
            Total.sd = sd(Total.coral),
            Total.SE = sd(Total.coral)/sqrt(n),
            BAcropora.mean = mean(Bleached.Acropora), 
            BAcropora.sd = sd(Bleached.Acropora),
            BAcropora.SE = sd(Bleached.Acropora)/sqrt(n()),
            BCoral.mean = mean(Bleached.other.coral), 
            BCoral.sd = sd(Bleached.other.coral),
            BCoral.SE = sd(Bleached.other.coral)/sqrt(n()),
            TotalB.mean = mean(Total.bleached), 
            TotalB.sd = sd(Total.bleached),
            TotalB.SE = sd(Total.bleached)/sqrt(n()),
            BTotal.mean = mean(Total.incl.bleached),
            BTotal.sd = sd(Total.incl.bleached),
            BTotal.SE = sd(Total.incl.bleached)/sqrt(n)) %>% 
    ungroup()
  
## SEs slightly different now as different sample size

## Mean by reef type ----

# From raw data:
mean_cover_rt_raw <- perc_all22 %>% 
  group_by(Reef.type) %>% 
  summarise(Acropora.mean = mean(Acropora), 
            Acropora.sd = sd(Acropora),
            n = n(),
            Acropora.SE = sd(Acropora)/sqrt(n()),
            Coral.mean = mean(Coral), 
            Coral.sd = sd(Coral),
            Coral.SE = sd(Coral)/sqrt(n()),
            Total.mean = mean(Total.coral),
            Total.sd = sd(Total.coral),
            Total.SE = sd(Total.coral)/sqrt(n)) %>%  
  ungroup()

# Try using function instead
rt_acro_mean <- summarySE(perc_all22, measurevar = "Acropora", groupvars = c("Reef.type"))

## gives same means, so much more efficient!!

# Plot mean Acropora cover by reef type

# Reorder categories
rt_acro_mean$Reef.type <- factor(rt_acro_mean$Reef.type,
                               levels = c("Channel", "Outer", "Inner", "Thila"),
                               labels = c("Channel", "Outer", "Inner", "Thila"))

# Plot as point and error bar
(rt_acro <- ggplot(rt_acro_mean, aes(x = Reef.type, y = Acropora)) +
  geom_point(size = 3, shape = 21, fill = "black") +
  geom_errorbar(aes(ymin = Acropora - se, ymax = Acropora + se), 
                width = .1) +
  expand_limits(y = 0) +
  theme_diss()
)

# Plot as bar plot
(rt_acro_bar <- ggplot(rt_acro_mean, aes(x = Reef.type, y = Acropora, 
                                         colour = Reef.type, fill = Reef.type)) +
    geom_bar(stat = "identity", colour = "black", size = .3, width = .7) +
    geom_errorbar(aes(ymax = Acropora + se, ymin = Acropora - se), 
                  width = .3, size = .3, colour = "black") +
    labs(y = "Mean Acropora cover (%)") +
    scale_y_continuous(limits = c(0,10),
                       expand = c(0,0)) + 
    scale_colour_manual(values = c('#225ea8','#41b6c4','#a1dab4','#edf8b1')) +
    scale_fill_manual(values = c('#225ea8','#41b6c4','#a1dab4','#edf8b1'))+
    theme_diss() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))
)

ggsave2(rt_acro_bar, filename = "outputs/rt_acro_bar.png", width = 5, height = 4)

## Mean by depth ----

dep_acro_mean <- summarySE(perc_all22, measurevar = "Acropora", groupvars = c("Depth"))

# Plot as bar plot
(dep_acro_bar <- ggplot(dep_acro_mean, aes(x = Depth, y = Acropora, 
                                         colour = Depth, fill = Depth)) +
    geom_bar(stat = "identity", colour = "black", size = .3, width = .6) +
    geom_errorbar(aes(ymax = Acropora + se, ymin = Acropora - se), 
                  width = .3, size = .3, colour = "black") +
    labs(y = "Mean Acropora cover (%)") +
    scale_y_continuous(limits = c(0,6),
                       expand = c(0,0)) + 
    scale_colour_manual(values = c('#2c7fb8','#7fcdbb','#edf8b1')) +
    scale_fill_manual(values = c('#2c7fb8','#7fcdbb','#edf8b1'))+
    theme_diss() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))
)

ggsave2(dep_acro_bar, filename = "outputs/dep_acro_bar.png", width = 5, height = 4)

## Mean by atoll side ----

EW_data <- perc_all22 %>% 
  filter(Reef.type %in% c("Channel", "Outer")) %>% 
  mutate(atoll.side = case_when(Site.ID %in% c("IM", "FK", "MBO", "MK", "FO") ~ "East",
                                Site.ID %in% c("VK", "MC", "MO", "FGO", "GO",
                                               "KO", "HC", "OC", "GC", "HW") ~ "West"))
# Reorder west before east
EW_data$atoll.side <- factor(EW_data$atoll.side, levels = c("West", "East"),
                                  labels = c("West", "East"))

ew_acro_mean <- summarySE(EW_data, measurevar = "Acropora", groupvars = c("atoll.side"))

# Plot as bar plot
(ew_acro_bar <- ggplot(ew_acro_mean, aes(x = atoll.side, y = Acropora, 
                                           colour = atoll.side, fill = atoll.side)) +
    geom_bar(stat = "identity", colour = "black", size = .3, width = .5) +
    geom_errorbar(aes(ymax = Acropora + se, ymin = Acropora - se), 
                  width = .2, size = .3, colour = "black") +
    labs(y = "Mean Acropora cover (%)") +
    scale_y_continuous(limits = c(0,8),
                       expand = c(0,0)) + 
    scale_colour_manual(values = c("#a1dab4", "#2c7fb8")) +
    scale_fill_manual(values = c("#a1dab4", "#2c7fb8"))+
    theme_diss() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))
)

ggsave2(ew_acro_bar, filename = "outputs/ew_acro_bar.png", width = 5, height = 4)


# Data exploration with data aggr. by depth ----

length(unique(all22sites$Site.name))  # 28 sites

# Mean Acropora cover: 
mean(all22sites$Acropora)  # 4.889018

# Mean other coral cover:
mean(all22sites$Other.Coral)  # 22.60464

# Mean total hard coral cover:
mean(all22sites$Total.hard.coral)  # 27.49366

# Percentage of acropora out of total:
4.889018/27.49366 * 100  # 17.78235 %

# Percentage of non-acropora out of total:
22.60464/27.49366 * 100  # 82.21765 %

### All the same averages, so will use full raw dataset for now ###


# Visualisation ----

## Stacked barplots ----

perc_long <- perc_all22 %>% 
  pivot_longer(cols = 7:25, names_to = "Cover", values_to = "Percent")

# With all categories
(stack <- ggplot(perc_long, aes(x = Depth, y = Percent, fill = Cover)) +
    geom_bar(position = "fill", stat = "identity", width = 0.4) +
    scale_y_continuous(labels = scales::percent,
                       expand = c(0,0)) + 
    #labs(y = "Benthic cover\n", x = "\nYear") +
    theme_diss() +
    theme(legend.position = "right",
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          #aspect.ratio = 1.5,
          panel.border = element_rect(fill = NA),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
)

# Filter out Acropora and other coral
coral_long <- perc_long %>% 
  filter(Cover == "Acropora" | Cover == "Coral")

# Add column to use as x-axis to show Laamu-wide average
coral_long <- coral_long %>% 
  mutate(X = "Laamu")

# Stacked barplot with only coral
(coral_stack <- ggplot(coral_long, aes(x = X, y = Percent, fill = Cover)) +
    geom_bar(position = "fill", stat = "identity", width = 0.5) +
    scale_y_continuous(labels = scales::percent,
                       expand = c(0,0)) +
    scale_colour_manual(values = c("#fcd0a1", "#b1b695")) +
    scale_fill_manual(values = c("#fcd0a1", "#b1b695")) +
    labs(y = "Percent cover") +
    theme_diss() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border = element_rect(fill = NA),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
)

save_plot(coral_stack, filename = "outputs/coral_stack.png")

ggsave(coral_stack, filename = "outputs/coral_stack.png", height = 6, width = 7)

# by depth
(coral_stack_dep <- ggplot(coral_long, aes(x = Depth, y = Percent, fill = Cover)) +
    geom_bar(position = "fill", stat = "identity", width = 0.4) +
    scale_y_continuous(labels = scales::percent,
                       expand = c(0,0)) + 
    scale_colour_manual(values = c("#fcd0a1", "#b1b695")) +
    scale_fill_manual(values = c("#fcd0a1", "#b1b695")) +
    labs(y = "Percent cover") +
    theme_diss() +
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          #axis.text.x = element_blank(),
          panel.border = element_rect(fill = NA),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
)
## very similar across depths

ggsave(coral_stack_dep, filename = "outputs/coral_stack_dep.png", height = 6, width = 7)

# Reorder reef types
coral_long$Reef.type <- factor(coral_long$Reef.type,
                                 levels = c("Channel", "Outer", "Inner", "Thila"),
                                 labels = c("Channel", "Outer", "Inner", "Thila"))

# by reef type
(coral_stack_rt <- ggplot(coral_long, aes(x = Reef.type, y = Percent, fill = Cover)) +
    geom_bar(position = "fill", stat = "identity", width = 0.4) +
    scale_y_continuous(labels = scales::percent,
                       expand = c(0,0)) + 
    scale_colour_manual(values = c("#fcd0a1", "#b1b695")) +
    scale_fill_manual(values = c("#fcd0a1", "#b1b695")) +
    labs(y = "Percent cover") +
    coord_flip() +
    theme_diss() +
    theme(legend.title = element_blank(),
          axis.title.y = element_blank(),
          #axis.text.x = element_blank(),
          panel.border = element_rect(fill = NA),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
)
# quite marked diffs, lowest in channel, and highest in thila
# why does it look higher for outer than inner reef here? Doesn't match with averages from earlier

# Plot with coord flipepd
(coral_stack_rt2 <- ggplot(coral_long, aes(x = Reef.type, y = Percent, fill = Cover)) +
    geom_bar(position = "fill", stat = "identity", width = 0.4) +
    scale_y_continuous(labels = scales::percent,
                       expand = c(0,0)) + 
    scale_colour_manual(values = c("#fcd0a1", "#b1b695")) +
    scale_fill_manual(values = c("#fcd0a1", "#b1b695")) +
    labs(y = "Percent cover") +
    coord_flip() +
    theme_diss() +
    theme(legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          #axis.text.x = element_blank(),
          panel.border = element_rect(fill = NA),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
)

# by site
(coral_stack_site <- ggplot(coral_long, aes(x = reorder(Site.name,(desc(Percent))), y = Percent, fill = Cover)) +
    geom_bar(position = "fill", stat = "identity", width = 0.4) +
    scale_y_continuous(labels = scales::percent,
                       expand = c(0,0)) + 
    labs(y = "Percent cover\n") +
    theme_diss() +
    theme(legend.position = "right",
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(fill = NA),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
)

# Arrange the data to have it ordered by cover
coral_long <- coral_long %>% 
  arrange(Cover != 'Acropora', Percent) %>%
  mutate(Site = factor(Site.name, unique(Site.name)),
         Cover = factor(Cover, unique(Cover)))
  
(coral_stack_site <- ggplot(coral_long, aes(x = Site, y = Percent, fill = Cover)) +
     geom_bar(position = "fill", stat = "identity", width = 0.4) +
     scale_y_continuous(labels = scales::percent,
                        expand = c(0,0)) + 
     labs(y = "Percent cover\n") +
     theme_diss() +
     theme(legend.position = "right",
           legend.title = element_blank(),
           axis.title.x = element_blank(),
           axis.text.x = element_text(angle = 45, hjust = 1),
           panel.border = element_rect(fill = NA),
           plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
)

coral_long %>%
  mutate(Site = fct_reorder(Site, Percent)) %>%
  ggplot(aes(x = Site, y = Percent, fill = Cover)) +
  geom_bar(position = "fill", stat = "identity", width = 0.4) +
  scale_y_continuous(labels = scales::percent,
                     expand = c(0,0)) +
  scale_colour_manual(values = c("#fcd0a1", "#b1b695")) +
  scale_fill_manual(values = c("#fcd0a1", "#b1b695")) +
  labs(y = "Percent cover") +
  theme_diss() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(fill = NA),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

coral_long %>%
  mutate(Site = fct_reorder(Site, Percent)) %>%
  ggplot(aes(x = Site, y = Percent, fill = Cover)) +
  geom_bar(position = "fill", stat = "identity", width = 0.4) +
  scale_y_continuous(labels = scales::percent,
                     expand = c(0,0)) +
  scale_colour_manual(values = c("#fcd0a1", "#b1b695")) +
  scale_fill_manual(values = c("#fcd0a1", "#b1b695")) +
  labs(y = "Percent cover") +
  coord_flip() +
  theme_diss() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(fill = NA),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))


## Pie chart of Laamu wide average ----

# Make subset with averaged proportions
coral_prop <- coral_long %>% 
  group_by(Cover) %>% 
  summarise(Prop = mean(Percent),
            Perc = Prop/27.49366) %>% 
  ungroup()

# With proportion
pie(coral_prop$Prop, labels = c("Acropora", "Other coral"))

# With percentage out of total
pie(coral_prop$Perc, labels = c("Acropora (17.8%)", "Other coral (82.2%)"), 
    col = c("#fcd0a1", "#b1b695"), border = "white")

## Visualise bleached out of total coral ----

# Remove total columns for visualisation processes
perc_rep_long <- perc_rep %>% 
  select(- Total.coral, - Total.bleached, - Total.incl.bleached)

perc_rep_long <- perc_rep_long %>% 
  pivot_longer(cols = 5:8, names_to = "Cover", values_to = "Percent")

# Add column to use as x-axis to show Laamu-wide average
perc_rep_long <- perc_rep_long %>% 
  mutate(X = "Laamu")

# Reorder categories
perc_rep_long$Cover <- factor(perc_rep_long$Cover,
                               levels = c("Acropora", "Bleached.Acropora",
                                          "Other.coral", "Bleached.other.coral"),
                               labels = c("Acropora", "Bleached.Acropora",
                                          "Other.coral", "Bleached.other.coral"))


(coral_stack <- ggplot(perc_rep_long, aes(x = X, y = Percent, fill = Cover)) +
    geom_bar(position = "fill", stat = "identity", width = 0.5) +
    scale_y_continuous(labels = scales::percent,
                       expand = c(0,0)) +
    scale_colour_manual(values = c("#cc8b86", "#f9eae1", "#7d4f50", "#d1be9c"),
                        labels = c("Acropora", "Bleached Acropora",
                                   "Other coral", "Bleached other coral")) +
    scale_fill_manual(values = c("#cc8b86", "#f9eae1", "#7d4f50", "#d1be9c"),
                      labels = c("Acropora", "Bleached Acropora",
                                 "Other coral", "Bleached other coral")) +
    labs(y = "Percent cover") +
    theme_diss() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.border = element_rect(fill = NA),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
)

# Same but as pie chart

# Make subset with averaged proportions
coral_prop <- perc_rep_long %>% 
  group_by(Cover) %>% 
  summarise(Prop = mean(Percent)) %>% 
  ungroup()

pie(coral_prop$Prop, labels = c("Acropora", "Bleached Acropora",
                                "Other coral", "Bleached other coral"),
    col = c("#cc8b86", "#f9eae1", "#7d4f50", "#d1be9c"), border = "white")

# Analysis ----

## Effect of depth on Acropora cover ----

hist(perc_rep$Acropora)

m1 <- lmer(Acropora ~ Depth + (1|Site.name), data = perc_rep)
summary(m1)
plot(m1)
simulateResiduals(m1, plot = T)  # red!
anova(m1)  # n.s

## Effect of reef type on Acropora cover ----

m2 <- lmer(Acropora ~ Reef.type + (1|Site.name), data = perc_rep)
summary(m2)
plot(m2)
simulateResiduals(m2, plot = T)  # red! (homogeneity of variance sig.)
anova(m2)  # n.s

## With raw data
m2.1 <- lmer(Acropora ~ Reef.type + (1|Site.name), data = perc_all22)
summary(m2.1)
plot(m2.1)
simulateResiduals(m2.1, plot = T) # very bad
anova(m2.1) n.s.

## Difference between sites ----

m3 <- lm(Acropora ~ Site.name, data = perc_rep)
summary(m3)
plot(m3)
simulateResiduals(m3, plot = T)  # bad fit

## Effect of atoll side on Acropora ----
m4 <- lmer(Acropora ~ atoll.side + (1|Site.name), data = EW_data)
summary(m4)
plot(m4)
simulateResiduals(m4, plot = T)  # very bad
anova(m4) n.s.
