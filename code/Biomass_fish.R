# 1.0 install packages ----
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("lme4")
install.packages("xlxs")
install.packages("pbkrtest")
install.packages("lmerTest")
install.packages("emmeans")

library("tidyverse")
library("ggplot2")
library("dplyr")
library("lme4")
library("xlxs")

## 1.1 - read in data ----
biom <- read.csv("fish_site_aggregated_data.csv")

# 2.0 - data tidying ----

biom1 <- biom %>%
  group_by(site, replicate_no, depth, species, family) %>%
  summarise(biom_sum = sum(biomass))

## 2.1 - export and re-read data for addition of reef type and importance ----
write.csv(biom1,"Desktop\\biom1.csv", row.names = TRUE)
biom1 <- read.csv("biom1_csv.csv")

# 3.0 numeric/factor/character changes ----
biom1$site <- as.factor(biom1$site)
biom1$reef_type <- as.factor(biom1$reef_type)
biom1$family <- as.factor(biom1$family)
biom1$species <- as.factor(biom1$species)
biom1$biom_sum <- as.numeric(biom1$biom_sum)
biom1$importance <- as.factor(biom1$importance)
biom1$biom_kg <- as.numeric(biom1$biom_kg)
biom1$depth <- as.factor(biom1$depth)

## data subsets
biom.gro <- filter(biom1, family %in% c("GRO"))
View(biom.gro)

--------------------------------------------------------------------------------------------------

# 4.0 - visualisations ----

## 4.1 - histogram of frequency of biomass ----
hist(biom1$biom_sum)

## 4.2 - com importance boxplot ----
com.biom <- filter(biom1, importance %in% c("com"))
View(com.biom)

com.plot<- ggplot(com.biom, aes(reef_type, biom_kg)) +
    geom_boxplot(fill = "#CD3333", alpha = 0.8, colour = "#8B2323") +  
    theme(axis.text.x = element_text(size = 12, angle = 0)) +
    labs(x = "reef_type", y = "biomass"))

## 4.2 - eco importance boxplot ----
eco.biom <- filter(biom1, importance %in% c("eco"))
View(eco.biom)

eco.plot<- ggplot(eco.biom, aes(reef_type, biom_kg)) +
    geom_boxplot(fill = "#CD3333", alpha = 0.8, colour = "#8B2323") +  
    theme(axis.text.x = element_text(size = 12, angle = 0)) +
    labs(x = "reef_type", y = "biomass"))

----------------------------------------------------------------------------------------------------------------------------------
# 5.0 - difference in biomass of ecologically important fish between reef type ----

## 5.1 - glm of biomass of eco important fish between reef type  ----

eco.glm <-glm(biom_sum ~ reef_type, data=eco.biom, family=quasipoisson(link="log"))
summary(eco.glm)

eco.null <- glm(biom_sum ~ 1, data=eco.biom, family=quasipoisson(link="log"))
summary(eco.null)

anova(eco.null, eco.glm, test="F")

## 5.2 - glm of biomass of com important fish between reef type  ----

com.glm <-glm(biom_sum ~ reef_type, data=com.biom, family=quasipoisson(link="log"))
summary(com.glm)

com.null <- glm(biom_sum ~ 1, data=com.biom, family=quasipoisson(link="log"))
summary(com.null)

anova(com.null, com.glm, test="F")

#model significantly improved by adding reef type

com.biom$reef_type = relevel(com.biom$reef_type, ref = "inside")

##higher for channel and outside reefs than inside reefs


## 5.3 - glm of biomass of groupers fish between reef type  ----

glm.gro <- glm(biom_kg ~ reef_type, data=gro.biom, family=quasipoisson(link="log"))
summary(glm.gro)
  
gro.null <- glm(biom_kg ~ 1, data=gro.biom, family=quasipoisson(link="log"))
summary(gro.null)

anova(gro.null, glm.gro, test="F")

## 5.4 - glm of biomass of snappers fish between reef type  ----

glm.snap <- glm(biom_kg ~ reef_type, data=snap.biom, family=quasipoisson(link="log"))
summary(glm.snap)

snap.null <- glm(biom_kg ~ 1, data=snap.biom, family=quasipoisson(link="log"))
summary(snap.null)

anova(snap.null, glm.snap, test="F")

## 5.5 - glm of biomass of emperors fish between reef type  ----

glm.emp <- glm(biom_kg ~ reef_type, data=emp.biom, family=quasipoisson(link="log"))
summary(glm.emp)

emp.null <- glm(biom_kg ~ 1, data=emp.biom, family=quasipoisson(link="log"))
summary(emp.null)

anova(emp.null, glm.emp, test="F")

## 5.6 - glm of biomass of parrotfish  between reef type  ----

par.biom <- filter(biom1, family %in% c("PAR"))
View(par.biom)

glm.par <- glm(biom_kg ~ reef_type, data=par.biom, family=quasipoisson(link="log"))
summary(glm.par)

par.null <- glm(biom_kg ~ 1, data=par.biom, family=quasipoisson(link="log"))
summary(par.null)

anova(par.null, glm.par, test="F")

## 5.7 - glm of biomass of rabbitfish  between reef type  ----

rab.biom <- filter(biom1, family %in% c("RAB"))
View(rab.biom)

glm.rab <- glm(biom_kg ~ reef_type, data=rab.biom, family=quasipoisson(link="log"))
summary(glm.rab)

rab.null <- glm(biom_kg ~ 1, data=rab.biom, family=quasipoisson(link="log"))
summary(rab.null)

anova(rab.null, glm.rab, test="F")

## 5.7 - glm of biomass of carangidae  between reef type  ----

jt.biom <- filter(biom1, family %in% c("T&J"))
View(jt.biom)

glm.jt <- glm(biom_kg ~ reef_type, data=jt.biom, family=quasipoisson(link="log"))
summary(glm.jt)

jt.null <- glm(biom_kg ~ 1, data=jt.biom, family=quasipoisson(link="log"))
summary(jt.null)

anova(jt.null, glm.jt, test="F")

jt.biom$reef_type = relevel(jt.biom$reef_type, ref = "inside")

----------------------------------------------------------------------------------------------------

### REDO WITH HELP FROM MEAGHAN ----

##notes from Meaghan: ----
#report emmeans by reporting by saying x has significantly more than y, 
#report estimate, z ratio and p value
#can report means and LCL and UCL in the sup material 

library(emmeans)  

# 1.0 Emperor linear mixed model ----

emp.lmm.reef <- lmer(biom_kg ~ reef_type + (1|site), data=emp.biom)
summary(emp.lmm.reef)  

## 1.1 emp null model
null.emp <- lmer(biom_kg ~ 1 + (1|site), data=emp.biom)
summary(null.emp)

anova(null.emp, emp.lmm.reef, test="F")

### p =  0.06097, df = 3, model not significant, chisq = 7.371

#emeans comparison between reef types for emperors

emmeans(emp.lmm.reef, pairwise ~ reef_type)

#no results are significant between reef types

---------------------------------------------------------------------
# 2.0 Snapper linear mixed model ----

snap.lmm <- lmer(biom_kg ~ reef_type + (1|site), data=snap.biom)
summary(snap.lmm)  

## 1.1 snap null model
null.snap <- lmer(biom_kg ~ 1 + (1|site), data=snap.biom)
summary(null.snap)

anova(null.snap, snap.lmm, test="F")

### p =  0.278, df = 3, model not significant, chisq = 3.851

#emeans comparison between reef types for emperors

emmeans(snap.lmm, pairwise ~ reef_type)

#no results are significant between reef types

---------------------------------------------------------------------
  # 3.0 Parrotfish linear mixed model ----

par.lmm <- lmer(biom_kg ~ reef_type + (1|site), data=par.biom)
summary(par.lmm)  

## 3.1 par null model
null.par <- lmer(biom_kg ~ 1 + (1|site), data=par.biom)
summary(null.par)

anova(null.par, par.lmm, test="F")

### p =  0.9124, df = 3, model not significant, chisq = 0.5295

#emeans comparison between reef types for emperors

emmeans(par.lmm, pairwise ~ reef_type)

#nothing is significant 

---------------------------------------------------------------------
  # 4.0 rabbitfish linear mixed model ----

rab.lmm <- lmer(biom_kg ~ reef_type + (1|site), data=rab.biom)
summary(rab.lmm)  

## 4.1 rab null model
null.rab <- lmer(biom_kg ~ 1 + (1|site), data=rab.biom)
summary(null.rab)

anova(null.rab, rab.lmm, test="F")

### p =  2.7405, df = 3, model not significant, chisq = 0.4334

#emeans comparison between reef types for emperors

emmeans(rab.lmm, pairwise ~ reef_type)

#nothing is significant 

---------------------------------------------------------------------
  # 5.0 T&J  linear mixed model ----

jt.lmm <- lmer(biom_kg ~ reef_type + (1|site), data=jt.biom)
summary(jt.lmm)  

## 5.1 t&j null model
null.jt <- lmer(biom_kg ~ 1 + (1|site), data=jt.biom)
summary(null.jt)

anova(null.jt, jt.lmm, test="F")

### p =  1.7698, df = 3, model not significant, chisq = 0.6215

#emeans comparison between reef types for emperors

emmeans(jt.lmm, pairwise ~ reef_type)

#nothing is significant 

---------------------------------------------------------------------
  # 6.0 grouper linear mixed model ----

gro.lmm <- lmer(biom_kg ~ reef_type + (1|site), data=gro.biom)
summary(gro.lmm)  

## 6.1 gro null model
null.gro <- lmer(biom_kg ~ 1 + (1|site), data=gro.biom)
summary(null.gro)

anova(null.gro, gro.lmm, test="F")

### p =  5.677, df = 3, model not significant, chisq = 0.1284

#emeans comparison between reef types for emperors

emmeans(gro.lmm, pairwise ~ reef_type)

#nothing is significant 

------------------------------------------------------------------
  
### NOTHING SIGNIFICANT BETWEEN REEF TYPES (CRY)
## TRYING BY DEPTH, SEE IF THAT WORKS 
  
------------------------------------------------------------------  
  
biom.lmm <- lmer(biom_kg ~ reef_type + (1|site), data=biom1)
summary(gro.lmm)  
  
  
  
  
  
  
  