# 1.0 install packages ----
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("lme4")
install.packages("pbkrtest")
install.packages("lmerTest")
install.packages("emmeans")

library("tidyverse")
library("ggplot2")
library("dplyr")
library("lme4")
library("emmeans")
---------------------------
# 2.1 certain columns ----
inverts <- invertebrate_data_sum_csv[,c(1, 2, 3, 4, 5, 6, 7)]
inverts <- inverts[-c(1473:1872),]

# 2.2 remove NAs ----
inverts <- inverts[is.na(inverts$site) == FALSE & is.na(inverts&spec) == FALSE
                   & is.na(inverts&year) == FALSE & is.na(inverts&fam) == FALSE
                   & is.na(inverts&depth) == FALSE,]

# 2.3 numeric/factor/character changes ----
inverts$site <- as.factor(inverts$site)
inverts$depth <- as.factor(inverts$depth)
inverts$rep <- as.factor(inverts$rep)
inverts$year <- as.factor(inverts$year)

----------------------------------------------------------
  
# 3.0 calculating abundance ----
inverts1 <- count(inverts, site, spec, depth, year, reef, rep)
inverts_year <- count(inverts, site, spec, year)  


# 4.1 boxplots abundance ----
(invert.plot <- ggplot(inverts1, aes(year, n)) +
   geom_boxplot(fill = "#CD3333", alpha = 0.8, colour = "#8B2323") +  
   theme(axis.text.x = element_text(size = 12, angle = 0)) +
   labs(x = "year", y = "n"))

# 4.2 checking for normality ----
hist(inverts1$n, main=NULL)
qqnorm(inverts1$n, main=NULL)
qqline(inverts1$n, main=NULL)
shapiro.test(inverts1$n)
#not normally distributed, W = 0.64292, p-value < 2.2e-16

-------------------------------------------------------------------------------------
  
# 5.1 - GLMM of species abundance by year with site as a random effect ----

invert.glmer <-glmer(n ~ year + (1|site), data=inverts1, family=poisson(link="log"))
summary(invert.glmer)
plot(invert.glmer)

#anova to compare to model that does not include year as an effect?
mod.null <- glmer(n ~ 1  + (1|site), data=inverts1, 
                  family=poisson(link="log"))
mod.null

anova(invert.glmer, mod.null, test="chisq")

----

# 5.2 - GLMM of species abundance by year including all other factors with site as a random effect ----
invert.glmer2 <-glmer(n ~ year + rep + depth + reef + (1|site), data=inverts1, family=poisson(link="log"))
summary(invert.glmer2)
plot(invert.glmer2)

#anova to compare to model that does not include year as an effect?
mod.null <- glmer(n ~ 1  + (1|site), data=inverts1, 
                  family=poisson(link="log"))
mod.null

anova(invert.glmer2, mod.null, test="chisq")


reduced.invert.1 <- update(invert.glmer2, . ~ . - year)
drop1(reduced.invert.1, test="Chi")

anova(reduced.invert.1, invert.glmer2, test="Chisq")
#reduced model explains slightly less variation but this is not statistically significant; can remove year
#from model without making it significantly worse


------------------------------------------------------------------------------------------------------

  # 6.0 - GLMM - effect of depth on abundance of invertebrates ----
  
inverts_depth <- count(inverts, site, spec, depth, year)
  
#### to do a boxplot of invert abundance by depth, facet wrapped by year 
#### run glmm of abundance by depth with site as random effect and year as random effect? 

invert.depth <-glmer(n ~ depth + year + (1|site), data=inverts_depth, family=poisson(link="log"))
summary(invert.depth)

null.invert.depth <- glmer(n ~ 1 + (1|site) , data=inverts_depth, family=poisson(link="log"))
summary(null.invert.depth)
anova(invert.depth, null.invert.depth, test="chisq")

#drop1 
reduced.depth.1 <- update(invert.depth, . ~ . - year)
drop1(reduced.depth.1, test="Chi")

emmeans(reduced.depth.1, pairwise ~ depth)


# 6.1 - boxplot depth as factor facet wrap by year

(invert.depth <- ggplot(inverts1, aes(depth, n)) +
    geom_boxplot(fill = "Light blue", alpha = 0.8, colour = "Dark blue") + facet_wrap("year") +
    theme(axis.text.x = element_text(size = 12, angle = 0)) +
    labs(x = "Depth", y = "Abundance")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



------------------------------------------------------------------------------------------------------
  
  # 7.0 - GLMM - effect of reef type on abundance of invertebrates ----
  
  inverts_reef <- count(inverts, site, spec, reef, year)

#### to do a boxplot of invert abundance by reef, facet wrapped by year 
#### run glmm of abundance by reef with site as random effect and year as random effect? 

invert.reef <-glmer(n ~ reef + year + (1|site), data=inverts_reef, family=poisson(link="log"))
summary(invert.reef)

#drop1 
reduced.reef.1 <- update(invert.reef, . ~ . - year)
drop1(reduced.reef.1, test="Chi")

#null mod
null.invert.reef <- glmer(n ~ 1 + (1|site) , data=inverts_reef, family=poisson(link="log"))
summary(null.invert.reef)

#anova
anova(reduced.reef.1, null.invert.reef, test="chisq")

#pairwise comparison
emmeans(reduced.reef.1, pairwise ~ reef)


# 7.1 - boxplot depth as factor facet wrap by year

(invert.reef.plot <- ggplot(inverts1, aes(reef, n)) +
    geom_boxplot(fill = "blue", alpha = 0.8, colour = "black") + 
    theme(axis.text.x = element_text(size = 12, angle = 0)) +
    labs(x = "reef", y = "n"))
  
  
------------------------------------------------------------------------
#8.0 breakdown by commercial/ecological importance----
inverts2 <- count(inverts, site, spec, depth, year, reef, rep, fam)

inverts2$importance[inverts2$fam == 'SC'] <- "com"
inverts2$importance[inverts2$fam == 'SEA'] <- "com"
inverts2$importance[inverts2$fam == 'LOBSTER'] <- "com"
inverts2$importance[inverts2$fam == 'GC'] <- "eco"
inverts2$importance[inverts2$fam == 'CUSHION'] <- "coral_feeder"
inverts2$importance[inverts2$fam == 'CS'] <- "coral_feeder"
inverts2$importance[inverts2$fam == 'COTS'] <- "coral_feeder"

----------------------------------------------------------------------------
  
# 9.0 - COM importance  ----
com.invert <- filter(inverts2, importance %in% c("com"))
View(com.invert)


## 9.1 boxplot commercially important by year ----

(year.com<- ggplot(com.invert, aes(year, n)) +
    geom_boxplot(fill = "blue", alpha = 0.8, colour = "black") + 
    theme(axis.text.x = element_text(size = 12, angle = 0)) +
    labs(x = "year", y = "n"))


## 9.2 run glmm of abundance by year with site as random effect and year as random effect ----

com.year.glmer <-glmer(n ~ depth + rep + year + reef + (1|site), data=com.invert, family=poisson(link="log"))
summary(com.year.glmer)

#drop1 
reduced.com.1 <- update(com.year.glmer, . ~ . - year)
drop1(reduced.com.1, test="Chi")

#null mod
null.com <- glmer(n ~ reef + rep + depth + (1|site) , data=com.invert, family=poisson(link="log"))
summary(null.com)

#anova
anova(com.year.glmer, null.com, test="chisq")


com.depth.glmer <-glmer(n ~ depth + rep + year + reef + (1|site), data=com.invert, family=poisson(link="log"))
summary(com.depth.glmer)

#null depth
null.depth.glmer <-glmer(n ~ rep + year + reef + (1|site), data=com.invert, family=poisson(link="log"))
summary(null.depth.glmer)

#anova
anova(com.depth.glmer, null.depth.glmer, test="chisq")

#pairwise comparison
emmeans(reduced.com.1, pairwise ~ year)
emmeans(com.depth.glmer, pairwise ~ depth)


----------------------------------------------------------------------------
  
  # 10.0 - coral feeding importance  ----
cf.invert <- filter(inverts2, importance %in% c("coral_feeder"))
View(cf.invert)


## 10.1 boxplot cf important by year ----

(year.cf<- ggplot(cf.invert, aes(year, n)) +
   geom_boxplot(fill = "blue", alpha = 0.8, colour = "black") + 
   theme(axis.text.x = element_text(size = 12, angle = 0)) +
   labs(x = "year", y = "n"))


## 10.2 run glmm of abundance by year with site as random effect and year as random effect ----

cf.year.glmer <-glmer(n ~ depth + rep + year + reef + (1|site), data=cf.invert, family=poisson(link="log"))
summary(cf.year.glmer)

#drop1 
reduced.cf.1 <- update(cf.year.glmer, . ~ . - reef)
drop1(reduced.cf.1, test="Chi")

#drop2
reduced.cf.2 <- update(reduced.cf.1, .~. - year)
drop1(reduced.cf.2, test="Chi")

#drop3
cf.year.glmer <-glmer(n ~ depth + rep + year + reef + (1|site), data=cf.invert, family=poisson(link="log"))
summary(cf.year.glmer)


#null mod
null.cf <- glmer(n ~ reef + depth + rep + (1|site) , data=cf.invert, family=poisson(link="log"))
summary(null.cf)

#anova
anova(cf.year.glmer, null.cf, test="chisq")
#model not improved by year p = 0.9289

-----
  
## 10.3 run glmm of abundance showing that reef type does not improve model

null.cf.reef <- glmer(n ~ year + depth + rep + (1|site) , data=cf.invert, family=poisson(link="log"))
summary(null.cf.reef)

#anova
anova(reduced.cf.1, null.cf.reef, test="chisq")
#model not improved by reef type p = 0.1818

-----
  
  ## 10.4 run glmm of abundance showing depth
  
  #
  cf.depth.glmer <-glmer(n ~ depth + (1|site), data=cf.invert, family=poisson(link="log"))
summary(cf.depth.glmer)


#null mod
null.cf <- glmer(n ~ 1 + (1|site) , data=cf.invert, family=poisson(link="log"))
summary(null.cf)

#anova
anova(cf.depth.glmer, null.cf, test="chisq")
#model significantly improved by depth p = 0.00249

#pairwise comparison
emmeans(cf.depth.glmer, pairwise ~ depth)



----------------------------------------------------------------------------
  
  # 11.0 - ecologically importance  ----
eco.invert <- filter(inverts2, importance %in% c("eco"))
View(eco.invert)


## 11.1 boxplot cf important by year ----

(year.eco<- ggplot(eco.invert, aes(year, n)) +
   geom_boxplot(fill = "blue", alpha = 0.8, colour = "black") + 
   theme(axis.text.x = element_text(size = 12, angle = 0)) +
   labs(x = "year", y = "n"))


## 11.2 run glmm of abundance by year with site as random effect and year as random effect ----

eco.year.glmer <-glmer(n ~ depth + year + reef + (1|site), data=eco.invert, family=poisson(link="log"))
summary(eco.year.glmer)

## 11.3 year null

null.eco.year <- glmer(n ~ + depth + reef + (1|site) , data=eco.invert, family=poisson(link="log"))
summary(null.eco.year)

#anova
anova(eco.year.glmer, null.eco.year, test="chisq")
#model improved by year p = 0.01743

----
  
## 11.4 run glmm of abundance by depth with site as random effect and year as random effect ----

eco.depth.glmer <-glmer(n ~ depth + year + reef + (1|site), data=eco.invert, family=poisson(link="log"))
summary(eco.depth.glmer)

## 11.5 depth null

null.eco.depth <- glmer(n ~ reef + year + (1|site) , data=eco.invert, family=poisson(link="log"))
summary(null.eco.depth)

#anova
anova(eco.depth.glmer, null.eco.depth, test="chisq")
#model improved by year p < 2.2e-16 ***

#pairwise comparison
emmeans(eco.depth.glmer, pairwise ~ depth)

## 11.6 reef analysis 

eco.reef.glmer <-glmer(n ~ depth + year + reef + (1|site), data=eco.invert, family=poisson(link="log"))
summary(eco.reef.glmer)

## 11.7 depth null

null.eco.reef <- glmer(n ~ depth + year + (1|site) , data=eco.invert, family=poisson(link="log"))
summary(null.eco.reef)
  
#anova
anova(eco.reef.glmer, null.eco.reef, test="chisq")
# model significantly improved by reef type p = 0.002483

#pairwise comparison
emmeans(eco.reef.glmer, pairwise ~ reef)


----------------------------------------------------------------------------
  
  # 12.0 - sea cucumbers importance  ----
inverts2$fam[inverts2$fam == 'SC'] <- "SEA"

sc.invert <- filter(inverts2, fam %in% c("SEA"))
View(sc.invert)

## 12.1 boxplot cf important by year ----

(year.sc<- ggplot(sc.invert, aes(year, n)) +
   geom_boxplot(fill = "blue", alpha = 0.8, colour = "black") + 
   theme(axis.text.x = element_text(size = 12, angle = 0)) +
   labs(x = "year", y = "n"))


## 12.2 run glmm of abundance by year with site as random effect and year as random effect ----

sc.year.glmer <-glmer(n ~ depth +  year + reef + (1|site), data=sc.invert, family=poisson(link="log"))
summary(sc.year.glmer)

#null no year 
null.sc.year <- glmer(n ~ depth + reef + (1|site) , data=sc.invert, family=poisson(link="log"))
summary(null.sc.year)

#anova
anova(sc.year.glmer, null.sc.year, test="chisq")

#model v significantly improved by year 

----
  ## 12.3 run glmm of abundance by year with site as random effect and year as random effect ----

sc.depth.glmer <-glmer(n ~ depth +  year + reef + (1|site), data=sc.invert, family=poisson(link="log"))
summary(sc.year.glmer)

#null no year 
null.sc.depth <- glmer(n ~ reef + year + (1|site) , data=sc.invert, family=poisson(link="log"))
summary(null.sc.year)

#anova
anova(sc.depth.glmer, null.sc.depth, test="chisq")

#model v significantly improved by year 

#pairwise comparison
emmeans(sc.depth.glmer, pairwise ~ depth)



