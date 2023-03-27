## Install Packages ----
#package installation ----
install.packages("ggplot2")
library("ggplot2")

install.packages("dplr")
library("dplyr")

install.packages("lme4")
library("lme4")


#run data and check it over ----
fam.abun <- family_abundance_year_site_com_eco[,c(1, 2, 3, 4, 5)]
summary(fam.abun)


#descriptive stats ----
table(fam.abun$site)
table(fam.abun$year)
table(fam.abun$abundance)
table(fam.abun$importance)
table(fam.abun$family)

counts.id <- table(fam.abun$importance)
counts.id

#remove NAs ----
fam.abun <- fam.abun[is.na(fam.abun$site)==
                         FALSE & is.na(fam.abun$family)==FALSE & is.na(fam.abun$abundance)
                       ==FALSE & is.na(fam.abun$site)
                     ==FALSE & is.na(fam.abun$year)==FALSE & is.na(fam.abun$importance)==FALSE,]

fam.abun

#numeric/factor/character changes ----
fam.abun$year <- as.factor(fam.abun$year)
fam.abun$site <- as.factor(fam.abun$site)
fam.abun$importance <- as.factor(fam.abun$importance)
fam.abun$abundance <- as.numeric(fam.abun$abundance)
fam.abun$family <- as.character(fam.abun$family)

#visualisations ----
##histogram of frequency of abundance ----
hist(fam.abun$abundance)


##table and barplot ----
table.fam.abun <- table(fam.abun$abundance)
barplot(table.fam.abun, xlab="number of fish families", ylab="frequency")

##boxplot ----
(fam.abun.plot <- ggplot(fam.abun, aes(year, abundance)) +
    geom_boxplot(fill = "#CD3333", alpha = 0.8, colour = "#8B2323") +  
    theme(axis.text.x = element_text(size = 12, angle = 0)) +
    labs(x = "year", y = "abundance"))

##data skewed by large value of gibbus at maavah_outside therefore removed from data and rerun above 

#Rerun without gibbus ----


#run data and check it over ----
fam.abun.nogib <- no_gib_family_abundance_year_site_com_eco[,c(1, 2, 3, 4, 5,6)]
summary(fam.abun.nogib)


#descriptive stats ----
table(fam.abun.nogib$site)
table(fam.abun.nogib$year)
table(fam.abun.nogib$abundance)
table(fam.abun.nogib$importance)
table(fam.abun.nogib$family)

counts.id <- table(fam.abun.nogib$importance)
counts.id

#remove NAs ----
fam.abun.nogib <- fam.abun.nogib[is.na(fam.abun.nogib$site)==
                       FALSE & is.na(fam.abun.nogib$family)==FALSE & is.na(fam.abun.nogib$abundance)
                     ==FALSE & is.na(fam.abun.nogib$site)
                     ==FALSE & is.na(fam.abun.nogib$year)==FALSE & is.na(fam.abun.nogib$importance)==FALSE,]

fam.abun.nogib

#numeric/factor/character changes ----
fam.abun.nogib$year <- as.factor(fam.abun.nogib$year)
fam.abun.nogib$site <- as.factor(fam.abun.nogib$site)
fam.abun.nogib$importance <- as.factor(fam.abun.nogib$importance)
fam.abun.nogib$abundance <- as.numeric(fam.abun.nogib$abundance)
fam.abun.nogib$family <- as.character(fam.abun.nogib$family)

#visualisations ----
##histogram of frequency of abundance ----
hist(fam.abun.nogib$abundance)


##table and barplot ----
table.fam.abun.nogib <- table(fam.abun.nogib$abundance)
barplot(table.fam.abun.nogib, xlab="number of fish families", ylab="frequency")

##boxplot ----
(fam.abun.plot.nogib <- ggplot(fam.abun.nogib, aes(year, abundance)) +
   geom_boxplot(fill = "#CD3333", alpha = 0.8, colour = "#8B2323") +  
   theme(axis.text.x = element_text(size = 12, angle = 0)) +
   labs(x = "year", y = "abundance"))

#checking for normality ----
par(mfrow=c(1,2))
hist(fam.abun.nogib$abundance, main=NULL)
qqnorm(fam.abun.nogib$abundance, main=NULL)
qqline(fam.abun.nogib$abundance, main=NULL)
shapiro.test(fam.abun.nogib$abundance)
#not normally distributed


#####################################################################################################################
# does year have an effect on the abundance by family? ----










 ----


##GLMM of family abundance by year with site as a random effect ############## ----

fam.glmer <-glmer(abundance ~ year + (1|site), data=fam.abun.nogib, family=poisson(link="log"))
summary(fam.glmer)
plot(fam.glmer)


##anova to compare to model that does not include year as an effect? ----
mod.null <- glmer(abundance ~ 1  + (1|site), data=fam.abun.nogib, 
                  family=poisson(link="log"))
mod.null

anova(mod.null, fam.glmer, test="chisq")

#significant difference in families between years


#####################################################################################################################
##does reef type have an effect on the abundance by family? ----


reef.fam.abun.glm.nogib <- glm(abundance ~ reef_type, data=fam.abun.nogib, 
                           family=quasipoisson(link="log"))
summary(reef.fam.abun.glm.nogib)

#looks as though reef type may effect the family abundance 

#anova to check
mod.null.reef.type <- glmer(abundance ~ 1  + (1|site), data=fam.abun.nogib, 
                  family=poisson(link="log"))
summary(mod.null.reef.type)

anova(mod.null.reef.type, reef.fam.abun.glm.nogib, test="chisq")


#####################################################################################################################
# abundance of commercially important vs ecologically important fish ----

##boxplot of commercially vs ecologically important fish for both years
dev.off()
(fam.abun.plot.nogib.importance <- ggplot(fam.abun.nogib, aes(importance, abundance)) +
     geom_boxplot(fill = "#CD3333", alpha = 0.8, colour = "#8B2323") +  
     theme(axis.text.x = element_text(size = 12, angle = 0)) +
     labs(x = "importance", y = "abundance"))

##com important fish abundance ----
levels(site.abun$year)
com.data <- filter(fam.abun.nogib, importance %in% c("com"))
View(com.data)

#boxplot of abundance of commercially important families
(com.year.plot.nogib.importance <- ggplot(com.data, aes(year, abundance)) +
        geom_boxplot(fill = "#CD3333", alpha = 0.8, colour = "#8B2323") +  
        theme(axis.text.x = element_text(size = 12, angle = 0)) +
        labs(x = "year", y = "abundance"))

#glmer of change in abundance of commercially important families between 2019-2022
fam.com.glmer <- glmer(abundance ~ year + (1|site), data=com.data, family=poisson(link="log"))
summary(fam.com.glmer)

plot(fam.com.glmer)

#anova 

mod.null.com <- glmer(abundance ~ 1  + (1|site), data=com.data, 
                            family=poisson(link="log"))
summary(mod.null.com)

anova(mod.null.com, fam.com.glmer, test="chisq")

##eco important fish abundance ----
levels(site.abun$year)
eco.data <- filter(fam.abun.nogib, importance %in% c("eco"))
View(eco.data)

#boxplot of abundance of commercially important families
(eco.year.plot.nogib.importance <- ggplot(eco.data, aes(year, abundance)) +
        geom_boxplot(fill = "#CD3333", alpha = 0.8, colour = "#8B2323") +  
        theme(axis.text.x = element_text(size = 12, angle = 0)) +
        labs(x = "year", y = "abundance"))


#glmer of change in abundance of ecologically important families between 2019-2022
fam.eco.glmer <- glmer(abundance ~ year + (1|site), data=eco.data, family=poisson(link="log"))
summary(fam.eco.glmer)

plot(fam.eco.glmer)

#anova 

mod.null.eco <- glmer(abundance ~ 1  + (1|site), data=eco.data, 
                      family=poisson(link="log"))
summary(mod.null.eco)

anova(fam.eco.glmer, mod.null.eco, test="chisq")

#####################################################################################################################
# effect of reef type on commercially important and ecologically important fish ----

##commercially important fish ----

##boxplot of commercially important fish for all reef types

(com.reef.plot.nogib <- ggplot(com.data, aes(reef_type, abundance)) +
     geom_boxplot(fill = "#CD3333", alpha = 0.8, colour = "#8B2323") +  
     theme(axis.text.x = element_text(size = 12, angle = 0)) +
     labs(x = "reef type", y = "abundance"))


#glm reef type, year and interaction with year
com.reef.year.abun.interaction.glm <- glm(abundance ~ reef_type + year + reef_type:year, data=com.data, 
                          family=quasipoisson(link="log"))
summary(com.reef.year.abun.interaction.glm)

#glm reef type and year
com.reef.type.year <- glm(abundance ~ reef_type + year, data=com.data, 
             family=quasipoisson(link="log"))
summary(com.reef.type.year)

#glm reef type alone
com.reef.glm <- glm(abundance ~ reef_type, data=com.data, 
                     family=quasipoisson(link="log"))
summary(com.reef.glm)


#null model
com.mod.null <- glm(abundance ~ 1, data=com.data, 
                family=quasipoisson(link="log"))
summary(com.mod.null)


#anova to compare null to models
##     
anova(com.mod.null,com.reef.year.abun.interaction.glm, test="F") #not significant

anova(com.mod.null, com.reef.type.year, test="F") #not significant

anova(com.mod.null, com.reef.glm, test="F") #not significant


##ecologically important fish ----

##boxplot of ecologically important fish for all reef types

(eco.reef.plot.nogib <- ggplot(eco.data, aes(reef_type, abundance)) +
     geom_boxplot(fill = "#CD3333", alpha = 0.8, colour = "#8B2323") +  
     theme(axis.text.x = element_text(size = 12, angle = 0)) +
     labs(x = "reef type", y = "abundance"))


#glm reef type, year and interaction with year
eco.reef.year.abun.interaction.glm <- glm(abundance ~ reef_type + year + reef_type:year, data=eco.data, 
                                          family=quasipoisson(link="log"))
summary(eco.reef.year.abun.interaction.glm)

#glm reef type and year
eco.reef.type.year <- glm(abundance ~ reef_type + year, data=eco.data, 
                          family=quasipoisson(link="log"))
summary(eco.reef.type.year)

#glm reef type alone
eco.reef.glm <- glm(abundance ~ reef_type, data=eco.data, 
                    family=quasipoisson(link="log"))
summary(eco.reef.glm)


#null model
eco.mod.null <- glm(abundance ~ 1, data=eco.data, 
                    family=quasipoisson(link="log"))
summary(eco.mod.null)


#anova to compare null to models
##     
anova(eco.mod.null,eco.reef.year.abun.interaction.glm, test="F") #not significant

anova(eco.mod.null, eco.reef.type.year, test="F") #not significant

anova(eco.mod.null, eco.reef.glm, test="F") #not significant






