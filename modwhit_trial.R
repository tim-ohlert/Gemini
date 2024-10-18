library(tidyverse)
library(plyr)
library(nlme)
library(emmeans)
library(vegan)
library(codyn)

modwhit <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/modwhit_clean_2024.csv")
Mod.whit.spp <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/Mod-whit-spp.csv")
transect.info <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/transect-id-2024.csv")
disturbance <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/modwhit-disturbance.csv")

all_but_thousand_m <- modwhit%>%
                  subset(Quad_sz_m2 != "1000" & Quad_sz_m2 != "1250")%>%
                  ddply(.(Year, Transect, Quad_num, Quad_sz_m2), function(x)data.frame(
                    species = ifelse(x$Spp_code == "none", 0,
                                length(x$Spp_code)
                  )))

thousand_m <- modwhit%>%
           subset(Spp_code != "none") %>%
ddply(.(Year, Transect), function(x)data.frame(
  species = length(unique(x$Spp_code))
))
thousand_m$Quad_num <- 14
thousand_m$Quad_sz_m2 <- ifelse(thousand_m$Transect == "B14" & thousand_m$Year == 2024, 1250, 1000) #one large plot got sampled at a larger area


four_thousand_m <- modwhit%>%
  left_join(transect.info, by = "Transect")%>%
  subset(Spp_code != "none") %>%
  ddply(.(Year, SoilVeg, Treatment), function(x)data.frame(
    species = length(unique(x$Spp_code))
  ))%>%
  subset(Year == 2024)

plot_summ <- rbind(all_but_thousand_m, thousand_m)%>%
              unique()%>%
              left_join(transect.info, by = "Transect")

plot_summ$log_sp <- log(plot_summ$species)
plot_summ$log_area <- log(plot_summ$Quad_sz_m2)

plot_summ <- dplyr::mutate(plot_summ, Treatment = fct_recode(Treatment, "Impact" = "Drive and Crush", "Control" = "Reference"), SoilVeg = fct_recode(SoilVeg, "Silty Saltbush" = "SiltyAtriplex", "Shallow Creosote" = "ShallowCreosote", "Deep Creosote" = "DeepCreosote"))



plot_summ%>%
  subset(Year == "2024")%>%
  ggplot(aes(log_area, log_sp, color = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()

plot_summ%>%
  subset(Year == "2024")%>%
  ggplot(aes(Quad_sz_m2, species, color = Treatment))+
  facet_wrap(~Transect)+
  geom_point()+
  geom_smooth(method = "glm", formula = y~log(x)#,
              #method.args = list(family = gaussian(link = 'log'))
              )+
  theme_bw()

plot_summ%>%
  subset(Year == "2024")%>%
  ggplot(aes(log_area, log_sp, color = Treatment))+
  facet_wrap(~Transect)+
  geom_point()+
  geom_smooth(method = "lm"#,
              #method.args = list(family = gaussian(link = 'log'))
  )+
  theme_bw()


    
mod <- lme(log_sp~log_area*Treatment, random = ~1|SoilVeg/Transect, data = plot_summ%>%
             subset(Year == "2024" & log_sp != "-Inf"))
summary(mod)


mod <- lm(species~log_area, data = plot_summ%>%
             subset(Transect == "B14" & Year == "2024" & log_sp != "-Inf"))
summary(mod)$coefficients[1]


metrics <- plot_summ%>%
  subset(Year == 2024 )%>%
  ddply(.(Transect, SoilVeg, Treatment), function(x)data.frame(
    intercept = summary(lm(log_sp~log_area, data = subset(x, log_sp != "-Inf")))$coefficients[1],
    slope = summary(lm(log_sp~log_area, data = subset(x, log_sp != "-Inf")))$coefficients[2]
  )) #add gamma diversity metric (saturation of curve)

#alpha diversity (intercept)

mod <- lme(intercept~Treatment, random = ~1|SoilVeg, data = metrics)
summary(mod)

mod <- lm(intercept~Treatment*SoilVeg, data = metrics)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))
#visreg(mod)

ggplot(metrics, aes(Treatment, intercept))+
  facet_wrap(~SoilVeg)+
  geom_point()+
  geom_boxplot()+
  ylab("Alpha diversity")+
  ylim(0,3)+
  theme_bw()

#beta diversity (slope)
mod <- lme(slope~Treatment, random = ~1|SoilVeg, data = metrics)
summary(mod)

mod <- lm(slope~Treatment*SoilVeg, data = metrics)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))
#visreg(mod)

ggplot(metrics, aes(Treatment, slope))+
  facet_wrap(~SoilVeg)+
  geom_point()+
  geom_boxplot()+
  ylab("Beta diversity")+
  theme_bw()


#gamma diversity (saturation point)


thousand_m_2023 <- rbind( thousand_m)%>%
  left_join(transect.info, by = "Transect")%>%
  dplyr::mutate( Treatment = fct_recode(Treatment, "Impact" = "Drive and Crush", "Control" = "Reference"), SoilVeg = fct_recode(SoilVeg, "Silty Saltbush" = "SiltyAtriplex", "Shallow Creosote" = "ShallowCreosote", "Deep Creosote" = "DeepCreosote"))


mod <- lme(species~Treatment, random = ~1|SoilVeg, data = thousand_m_2023)
summary(mod)


mod <- lm(species~Treatment*SoilVeg, data = thousand_m_2023)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))
#visreg(mod)


ggplot(thousand_m_2023, aes(Treatment, species))+
  facet_wrap(~SoilVeg)+
  geom_point()+
  geom_boxplot()+
  ylab("Gamma diversity")+
  ylim(0, 50)+
  theme_bw()


ggplot(four_thousand_m, aes(Treatment, species))+
  facet_wrap(~SoilVeg)+
  geom_point()+
  ylab("4km squared diversity")+
  ylim(0, 60)+
  theme_bw()


###### native species

all_but_thousand_m <- modwhit%>%
  left_join(Mod.whit.spp, by = join_by("Spp_code" == "SPP"))%>%
  subset(native == "Native")%>%
  subset(Quad_sz_m2 != "1000" & Quad_sz_m2 != "1250")%>%
  ddply(.(Year, Transect, Quad_num, Quad_sz_m2), function(x)data.frame(
    species = ifelse(x$Spp_code == "none", 0,
                     length(x$Spp_code)
    )))

thousand_m <- modwhit%>%
  left_join(Mod.whit.spp, by = join_by("Spp_code" == "SPP"))%>%
  subset(native == "Native")%>%
  subset(Spp_code != "none") %>%
  ddply(.(Year, Transect), function(x)data.frame(
    species = length(unique(x$Spp_code))
  ))
thousand_m$Quad_num <- 14
thousand_m$Quad_sz_m2 <- ifelse(thousand_m$Transect == "B14" & thousand_m$Year == 2024, 1250, 1000) #one large plot got sampled at a larger area




plot_summ <- rbind(all_but_thousand_m, thousand_m)%>%
  unique()%>%
  left_join(transect.info, by = "Transect")

plot_summ$log_sp <- log(plot_summ$species)
plot_summ$log_area <- log(plot_summ$Quad_sz_m2)

#plot_summ <- dplyr::mutate(plot_summ, Treatment = fct_recode(Treatment, "Impact" = "Drive and Crush", "Control" = "Reference"), SoilVeg = fct_recode(SoilVeg, "Silty Saltbush" = "SiltyAtriplex", "Shallow Creosote" = "ShallowCreosote", "Deep Creosote" = "DeepCreosote"))



plot_summ%>%
  subset(Year == "2024")%>%
  ggplot(aes(log_area, log_sp, color = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()

plot_summ%>%
  subset(Year == "2024")%>%
  ggplot(aes(Quad_sz_m2, species, color = Treatment))+
  facet_wrap(~Transect)+
  geom_point()+
  geom_smooth(method = "glm", formula = y~log(x)#,
              #method.args = list(family = gaussian(link = 'log'))
  )+
  theme_bw()

plot_summ%>%
  subset(Year == "2024")%>%
  ggplot(aes(log_area, log_sp, color = Treatment))+
  facet_wrap(~Transect)+
  geom_point()+
  geom_smooth(method = "lm"#,
              #method.args = list(family = gaussian(link = 'log'))
  )+
  theme_bw()



#mod <- lme(log_sp~log_area*Treatment, random = ~1|SoilVeg/Transect, data = plot_summ%>%
#             subset(Year == "2024" & log_sp != "-Inf"))
#summary(mod)


#mod <- lm(species~log_area, data = plot_summ%>%
#            subset(Transect == "B14" & Year == "2024" & log_sp != "-Inf"))
#summary(mod)$coefficients[1]


metrics <- plot_summ%>%
  subset(Year == 2024 )%>%
  ddply(.(Transect, SoilVeg, Treatment), function(x)data.frame(
    intercept = summary(lm(log_sp~log_area, data = subset(x, log_sp != "-Inf")))$coefficients[1],
    slope = summary(lm(log_sp~log_area, data = subset(x, log_sp != "-Inf")))$coefficients[2]
  )) #add gamma diversity metric (saturation of curve)

#alpha diversity (intercept)

mod <- lme(intercept~Treatment, random = ~1|SoilVeg, data = metrics)
summary(mod)

mod <- lm(intercept~Treatment*SoilVeg, data = metrics)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))


ggplot(metrics, aes(Treatment, intercept))+
  facet_wrap(~SoilVeg)+
  geom_point()+
  geom_boxplot()+
  ylab("Alpha diversity")+
  ylim(0,3)+
  theme_bw()

#beta diversity (slope)
mod <- lme(slope~Treatment, random = ~1|SoilVeg, data = metrics)
summary(mod)

mod <- lm(slope~Treatment*SoilVeg, data = metrics)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))
#visreg(mod)

ggplot(metrics, aes(Treatment, slope))+
  facet_wrap(~SoilVeg)+
  geom_point()+
  geom_boxplot()+
  ylab("Beta diversity")+
  theme_bw()


#gamma diversity (saturation point)


thousand_m_2023 <- rbind( thousand_m)%>%
  left_join(transect.info, by = "Transect")%>%
  dplyr::mutate( Treatment = fct_recode(Treatment, "Impact" = "Drive and Crush", "Control" = "Reference"), SoilVeg = fct_recode(SoilVeg, "Silty Saltbush" = "SiltyAtriplex", "Shallow Creosote" = "ShallowCreosote", "Deep Creosote" = "DeepCreosote"))


mod <- lme(species~Treatment, random = ~1|SoilVeg, data = thousand_m_2023)
summary(mod)


mod <- lm(species~Treatment*SoilVeg, data = thousand_m_2023)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))
#visreg(mod)


ggplot(thousand_m_2023, aes(Treatment, species))+
  facet_wrap(~SoilVeg)+
  geom_point()+
  geom_boxplot()+
  ylab("Gamma diversity")+
  ylim(0, 50)+
  theme_bw()



########### Annual species only (i.e. species most likely to have responded to disturbance)

all_but_thousand_m <- modwhit%>%
  left_join(Mod.whit.spp, by = join_by("Spp_code" == "SPP"))%>%
  subset(growth == "Annual")%>%
  subset(Quad_sz_m2 != "1000" & Quad_sz_m2 != "1250")%>%
  ddply(.(Year, Transect, Quad_num, Quad_sz_m2), function(x)data.frame(
    species = ifelse(x$Spp_code == "none", 0,
                     length(x$Spp_code)
    )))

thousand_m <- modwhit%>%
  left_join(Mod.whit.spp, by = join_by("Spp_code" == "SPP"))%>%
  subset(growth == "Annual")%>%
  subset(Spp_code != "none") %>%
  ddply(.(Year, Transect), function(x)data.frame(
    species = length(unique(x$Spp_code))
  ))
thousand_m$Quad_num <- 14
thousand_m$Quad_sz_m2 <- ifelse(thousand_m$Transect == "B14" & thousand_m$Year == 2024, 1250, 1000) #one large plot got sampled at a larger area




plot_summ <- rbind(all_but_thousand_m, thousand_m)%>%
  unique()%>%
  left_join(transect.info, by = "Transect")

plot_summ$log_sp <- log(plot_summ$species)
plot_summ$log_area <- log(plot_summ$Quad_sz_m2)



plot_summ%>%
  subset(Year == "2024")%>%
  ggplot(aes(log_area, log_sp, color = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()

plot_summ%>%
  subset(Year == "2024")%>%
  ggplot(aes(Quad_sz_m2, species, color = Treatment))+
  facet_wrap(~Transect)+
  geom_point()+
  geom_smooth(method = "glm", formula = y~log(x)#,
              #method.args = list(family = gaussian(link = 'log'))
  )+
  theme_bw()

plot_summ%>%
  subset(Year == "2024")%>%
  ggplot(aes(log_area, log_sp, color = Treatment))+
  facet_wrap(~Transect)+
  geom_point()+
  geom_smooth(method = "lm"#,
              #method.args = list(family = gaussian(link = 'log'))
  )+
  theme_bw()



metrics <- plot_summ%>%
  subset(Year == 2024 )%>%
  ddply(.(Transect, SoilVeg, Treatment), function(x)data.frame(
    intercept = summary(lm(log_sp~log_area, data = subset(x, log_sp != "-Inf")))$coefficients[1],
    slope = summary(lm(log_sp~log_area, data = subset(x, log_sp != "-Inf")))$coefficients[2]
  )) #add gamma diversity metric (saturation of curve)

#alpha diversity (intercept)

mod <- lme(intercept~Treatment, random = ~1|SoilVeg, data = metrics)
summary(mod)

mod <- lm(intercept~Treatment*SoilVeg, data = metrics)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))


ggplot(metrics, aes(Treatment, intercept))+
  facet_wrap(~SoilVeg)+
  geom_point()+
  geom_boxplot()+
  ylab("Alpha diversity")+
  ylim(0,3)+
  theme_bw()

#beta diversity (slope)
mod <- lme(slope~Treatment, random = ~1|SoilVeg, data = metrics)
summary(mod)

mod <- lm(slope~Treatment*SoilVeg, data = metrics)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))
#visreg(mod)

ggplot(metrics, aes(Treatment, slope))+
  facet_wrap(~SoilVeg)+
  geom_point()+
  geom_boxplot()+
  ylab("Beta diversity")+
  theme_bw()


#gamma diversity (saturation point)


thousand_m_2023 <- rbind( thousand_m)%>%
  left_join(transect.info, by = "Transect")%>%
  dplyr::mutate( Treatment = fct_recode(Treatment, "Impact" = "Drive and Crush", "Control" = "Reference"), SoilVeg = fct_recode(SoilVeg, "Silty Saltbush" = "SiltyAtriplex", "Shallow Creosote" = "ShallowCreosote", "Deep Creosote" = "DeepCreosote"))


mod <- lme(species~Treatment, random = ~1|SoilVeg, data = thousand_m_2023)
summary(mod)


mod <- lm(species~Treatment*SoilVeg, data = thousand_m_2023)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))
#visreg(mod)


ggplot(thousand_m_2023, aes(Treatment, species))+
  facet_wrap(~SoilVeg)+
  geom_point()+
  geom_boxplot()+
  ylab("Gamma diversity")+
  ylim(0, 50)+
  theme_bw()




#######try out metacommunity difference analysis

#all_but_thousand_m <- modwhit%>%
#  left_join(Mod.whit.spp, by = join_by("Spp_code" == "SPP"))%>%
#  subset(growth == "Annual")%>%
#  subset(Quad_sz_m2 != "1000" & Quad_sz_m2 != "1250")%>%
#  ddply(.(Year, Transect, Quad_num, Quad_sz_m2), function(x)data.frame(
#    species = ifelse(x$Spp_code == "none", 0,
#                     length(x$Spp_code)
#    )))

thousand_m <- modwhit%>%
  left_join(Mod.whit.spp, by = join_by("Spp_code" == "SPP"))%>%
  left_join(transect.info, by = "Transect")%>%
  #subset(growth == "Annual")%>%
  subset(Spp_code != "none" ) %>%
  subset(Year == 2024)%>%
  dplyr::select(Transect, Spp_code, SoilVeg, Treatment)%>%
  unique()
  #ddply(.(Year, Transect), function(x)data.frame(
  #  species = length(unique(x$Spp_code))
 # ))
#thousand_m$Quad_num <- 14
#thousand_m$Quad_sz_m2 <- ifelse(thousand_m$Transect == "B14" & thousand_m$Year == 2024, 1250, 1000) #one large plot got sampled at a larger area


#spread
thousand_m$present <- 1
thousand.spread <- pivot_wider(thousand_m, names_from = Spp_code, values_from = present, values_fill = 0)

library(vegan)
#
matrix_1 <- thousand.spread[,7:86]
matrix <- sapply( matrix_1, as.numeric )
## Bray-Curtis distances between samples
dis <- vegdist(matrix)
## extract treatments
groups <- paste(thousand.spread$Treatment, thousand.spread$SoilVeg)
## Calculate multivariate dispersions
mod <- betadisper(dis, groups, type = "centroid")
mod
## Perform test
anova(mod)
## Permutation test for F
permutest(mod, pairwise = TRUE, permutations = 99)

## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(mod)
## with data ellipses instead of hulls
plot(mod, ellipse = TRUE, hull = FALSE) # 1 sd data ellipse
plot(mod, ellipse = TRUE, hull = FALSE, conf = 0.90) # 90% data ellipse
## can also specify which axes to plot, ordering respected
plot(mod, axes = c(1,2), seg.col = "forestgreen", seg.lty = "dashed")
## Draw a boxplot of the distances to centroid for each group
boxplot(mod)



######2023 Deep creosote check

thousand_m <- modwhit%>%
  left_join(Mod.whit.spp, by = join_by("Spp_code" == "SPP"))%>%
  left_join(transect.info, by = "Transect")%>%
  #subset(growth == "Annual")%>%
  subset(Spp_code != "none" ) %>%
  subset(Year == 2023)%>%
  dplyr::select(Transect, Spp_code, SoilVeg, Treatment)%>%
  unique()

#spread
thousand_m$present <- 1
thousand.spread <- pivot_wider(thousand_m, names_from = Spp_code, values_from = present, values_fill = 0)

library(vegan)
#
matrix_1 <- thousand.spread[,7:66]
matrix <- sapply( matrix_1, as.numeric )
## Bray-Curtis distances between samples
dis <- vegdist(matrix)
## extract treatments
groups <- paste(thousand.spread$Treatment, thousand.spread$SoilVeg)
## Calculate multivariate dispersions
mod <- betadisper(dis, groups, type = "centroid")
mod
## Perform test
anova(mod)
## Permutation test for F
permutest(mod, pairwise = TRUE, permutations = 99)

## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(mod)
## with data ellipses instead of hulls
plot(mod, ellipse = TRUE, hull = FALSE) # 1 sd data ellipse
plot(mod, ellipse = TRUE, hull = FALSE, conf = 0.90) # 90% data ellipse
## can also specify which axes to plot, ordering respected
plot(mod, axes = c(1,2), seg.col = "forestgreen", seg.lty = "dashed")
## Draw a boxplot of the distances to centroid for each group
boxplot(mod)



####rank abundance
rank_abun <- modwhit%>%
  subset(Quad_sz_m2 == 1 & Spp_code != "none")%>%
  subset(Year == 2024)%>%
  dplyr::select(Transect, Quad_num,Spp_code, Cover_perc)%>%
  pivot_wider(names_from = "Spp_code", values_from = "Cover_perc", values_fill = 0)%>%
  pivot_longer(3:62, names_to = "species", values_to = "cover")%>%
  left_join(transect.info, by = "Transect")%>%
  group_by(SoilVeg, Treatment, species)%>%
  dplyr::summarize(se = sd(cover,na.rm = TRUE)/sqrt(length(cover)),cover = mean(cover, na.rm = TRUE) )

rank_abun%>%
  #subset(cover > 0.5)%>%
  subset(species == "AMDU2" | species == "ANLA7" |species == "CHRI" |species == "KRER" |species == "LATR2" |species == "MAAF" |species == "SCHIS")%>%
ggplot( aes(species, cover))+
  facet_grid(SoilVeg~Treatment)+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(x = species, ymin= cover-se, ymax = cover+se))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))



## rank abundance change
rank_abun_change <- modwhit%>%
  subset(Quad_sz_m2 == 1 & Spp_code != "none")%>%
  dplyr::select(Year, Transect, Quad_num,Spp_code, Cover_perc)%>%
  pivot_wider(names_from = "Spp_code", values_from = "Cover_perc", values_fill = 0)%>%
  pivot_longer(4:68, names_to = "species", values_to = "cover")%>%
  left_join(transect.info, by = "Transect")%>%
  pivot_wider(names_from = "Year", values_from = "cover")

rank_abun_change$change <- rank_abun_change$'2024' - rank_abun_change$'2023'

rank_abun_change <- rank_abun_change%>%
  group_by( SoilVeg, Treatment, Transect, species)%>%
  dplyr::summarize(change = mean(change, na.rm = TRUE) )%>%
  group_by( SoilVeg, Treatment, species)%>%
  dplyr::summarize(se = sd(change,na.rm = TRUE)/sqrt(length(change)),cover = mean(change, na.rm = TRUE) )

rank_abun_change_diff <- rank_abun_change%>%
  dplyr::select(SoilVeg, Treatment, species, cover)%>%
  pivot_wider(names_from = "Treatment", values_from = "cover")

rank_abun_change_diff$diff <- rank_abun_change_diff$'Drive and Crush' - rank_abun_change_diff$Reference


rank_abun_change%>%
  subset(cover > 0.5 | cover <0.5)%>%
  #subset(species == "AMDU2" | species == "ANLA7" |species == "CHRI" |species == "KRER" |species == "LATR2" |species == "MAAF" |species == "SCHIS")%>%
  ggplot( aes(species, cover))+
  facet_grid(SoilVeg~Treatment)+
  #geom_point()+
  geom_bar(stat = "identity")+
  ylab("Change in cover")+
  geom_errorbar(aes(x = species, ymin= cover-se, ymax = cover+se))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))




rank_abun_change_diff%>%
  left_join(Mod.whit.spp, join_by(species == SPP))%>%
  subset(diff > 0.25 | diff < -0.25)%>%
  subset(SoilVeg == "DeepCreosote")%>%
  #subset(species == "AMDU2" | species == "ANLA7" |species == "CHRI" |species == "KRER" |species == "LATR2" |species == "MAAF" |species == "SCHIS")%>%
  ggplot( aes(taxa, diff))+
#  facet_wrap(~SoilVeg)+
  #geom_point()+
  geom_bar(stat = "identity")+
  ylab("Change in cover relative to reference")+
  geom_hline(yintercept = 0)+
  #geom_errorbar(aes(x = species, ymin= diff-se, ymax = diff+se))+
  ylim(-3.5, 3.6)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))


rank_abun_change_diff%>%
  left_join(Mod.whit.spp, join_by(species == SPP))%>%
  subset(diff > 0.25 | diff < -0.25)%>%
  subset(SoilVeg == "ShallowCreosote")%>%
  #subset(species == "AMDU2" | species == "ANLA7" |species == "CHRI" |species == "KRER" |species == "LATR2" |species == "MAAF" |species == "SCHIS")%>%
  ggplot( aes(taxa, diff))+
  #  facet_wrap(~SoilVeg)+
  #geom_point()+
  geom_bar(stat = "identity")+
  ylab("Change in cover relative to reference")+
  geom_hline(yintercept = 0)+
  #geom_errorbar(aes(x = species, ymin= diff-se, ymax = diff+se))+
  ylim(-3.5, 3.6)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))





rank_abun_change_diff%>%
  left_join(Mod.whit.spp, join_by(species == SPP))%>%
  subset(diff > 0.25 | diff < -0.25)%>%
  subset(SoilVeg == "SiltyAtriplex")%>%
  #subset(species == "AMDU2" | species == "ANLA7" |species == "CHRI" |species == "KRER" |species == "LATR2" |species == "MAAF" |species == "SCHIS")%>%
  ggplot( aes(taxa, diff))+
  #  facet_wrap(~SoilVeg)+
  #geom_point()+
  geom_bar(stat = "identity")+
  ylab("Change in cover relative to reference")+
  geom_hline(yintercept = 0)+
  #geom_errorbar(aes(x = species, ymin= diff-se, ymax = diff+se))+
  ylim(-3.5, 3.6)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))



#############
###species changes across disturbance gradients
diff <- modwhit%>%
  subset(Quad_sz_m2 == 1 & Spp_code != "none")%>%
  dplyr::select(Year, Transect, Quad_num,Spp_code, Cover_perc)%>%
  pivot_wider(names_from = "Spp_code", values_from = "Cover_perc", values_fill = 0)%>%
  pivot_longer(4:68, names_to = "species", values_to = "cover")%>%
  pivot_wider(names_from = "Year", values_from = "cover")
  #left_join(transect.info, by = "Transect")

diff$change <- diff$'2023' - diff$'2024'

diff <- left_join(diff, disturbance, by = "Transect")

diff%>%
  subset(Treatment != "Control")%>%
ggplot( aes(perc_disturbance_2024, change, color = SoilVeg))+
  facet_wrap(~species)+
  geom_point()+
  geom_smooth()+
  ylim(-25, 25)+
  theme_bw()




############
################
##compare alpha, beta, gamma over time

metrics <- plot_summ%>%
  ddply(.(Year,Transect, SoilVeg, Treatment), function(x)data.frame(
    intercept = summary(lm(log_sp~log_area, data = subset(x, log_sp != "-Inf")))$coefficients[1],
    slope = summary(lm(log_sp~log_area, data = subset(x, log_sp != "-Inf")))$coefficients[2]
  )) #add gamma diversity metric (saturation of curve)

#alpha diversity (intercept)

mod <- lme(intercept~Treatment*as.factor(Year), random = ~1|SoilVeg, data = metrics)
summary(mod)

mod <- lm(intercept~Treatment*SoilVeg*as.factor(Year), data = metrics)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))
#visreg(mod)

ggplot(metrics, aes(Treatment, intercept, color = as.factor(Year)))+
  facet_wrap(~SoilVeg)+
  geom_point()+
  geom_boxplot()+
  ylab("Alpha diversity")+
  ylim(0,3)+
  theme_bw()

#beta diversity (slope)
mod <- lme(slope~Treatment*as.factor(Year), random = ~1|SoilVeg, data = metrics)
summary(mod)

mod <- lm(slope~Treatment*SoilVeg*as.factor(Year), data = metrics)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))
#visreg(mod)

ggplot(metrics, aes(Treatment, slope, color = as.factor(Year)))+
  facet_wrap(~SoilVeg)+
  geom_point()+
  geom_boxplot()+
  ylab("Beta diversity")+
  theme_bw()


#gamma diversity
thousand_m <- modwhit%>%
  subset(Spp_code != "none") %>%
  ddply(.(Year, Transect), function(x)data.frame(
    species = length(unique(x$Spp_code))
  ))
thousand_m$Quad_num <- 14
thousand_m$Quad_sz_m2 <- ifelse(thousand_m$Transect == "B14" & thousand_m$Year == 2024, 1250, 1000) #one large plot got sampled at a larger area
thousand_m_all <- rbind( thousand_m)%>%
  left_join(transect.info, by = "Transect")#%>%
#  dplyr::mutate( Treatment = fct_recode(Treatment, "Impact" = "Drive and Crush", "Control" = "Reference"), SoilVeg = fct_recode(SoilVeg, "Silty Saltbush" = "SiltyAtriplex", "Shallow Creosote" = "ShallowCreosote", "Deep Creosote" = "DeepCreosote"))


mod <- lme(species~Treatment*as.factor(Year), random = ~1|SoilVeg, data = thousand_m_all)
summary(mod)


mod <- lm(species~Treatment*SoilVeg*as.factor(Year), data = thousand_m_all)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))
#visreg(mod)


ggplot(thousand_m_all, aes(Treatment, species, color = as.factor(Year)))+
  facet_wrap(~SoilVeg)+
  geom_point()+
  geom_boxplot()+
  ylab("Gamma diversity")+
  ylim(0, 50)+
  theme_bw()


ggplot(four_thousand_m, aes(Treatment, species, color = as.factor(Year)))+
  facet_wrap(~SoilVeg)+
  geom_point()+
  ylab("4km squared diversity")+
  ylim(0, 60)+
  theme_bw()


#################
########Differences between 2023-2024

all_but_thousand_m <- modwhit%>%
  subset(Quad_sz_m2 != "1000" & Quad_sz_m2 != "1250")%>%
  ddply(.(Year, Transect, Quad_num, Quad_sz_m2), function(x)data.frame(
    species = ifelse(x$Spp_code == "none", 0,
                     length(x$Spp_code)
    )))

thousand_m <- modwhit%>%
  subset(Spp_code != "none") %>%
  ddply(.(Year, Transect), function(x)data.frame(
    species = length(unique(x$Spp_code))
  ))
thousand_m$Quad_num <- 14
thousand_m$Quad_sz_m2 <- ifelse(thousand_m$Transect == "B14" & thousand_m$Year == 2024, 1250, 1000) #one large plot got sampled at a larger area


#four_thousand_m <- modwhit%>%
#  left_join(transect.info, by = "Transect")%>%
#  subset(Spp_code != "none") %>%
#  ddply(.(Year, SoilVeg, Treatment), function(x)data.frame(
#    species = length(unique(x$Spp_code))
#  ))



plot_summ <- rbind(all_but_thousand_m, thousand_m)%>%
  unique()%>%
  left_join(transect.info, by = "Transect")

plot_summ$log_sp <- log(plot_summ$species)
plot_summ$log_area <- log(plot_summ$Quad_sz_m2)

plot_summ <- dplyr::mutate(plot_summ, Treatment = fct_recode(Treatment, "Impact" = "Drive and Crush", "Control" = "Reference"), SoilVeg = fct_recode(SoilVeg, "Silty Saltbush" = "SiltyAtriplex", "Shallow Creosote" = "ShallowCreosote", "Deep Creosote" = "DeepCreosote"))

diff_metrics <- plot_summ%>%
  #subset(Year == 2024 )%>%
  ddply(.(Year, Transect, SoilVeg, Treatment), function(x)data.frame(
    intercept = summary(lm(log_sp~log_area, data = subset(x, log_sp != "-Inf")))$coefficients[1],
    slope = summary(lm(log_sp~log_area, data = subset(x, log_sp != "-Inf")))$coefficients[2]
  ))%>%
  pivot_wider(names_from = Year, values_from = c("intercept", "slope"))%>%
  mutate(intercept_diff = intercept_2024-intercept_2023, slope_diff = slope_2024-slope_2023)

##alpha diversity
ggplot(diff_metrics, aes(Treatment, intercept_diff))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  geom_hline(yintercept = 0)+
  ylab("Change in alpha diversity (intercept)")+
  theme_bw()

##beta diversity
ggplot(diff_metrics, aes(Treatment, slope_diff))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  geom_hline(yintercept = 0)+
  ylab("Change in beta diversity (slope)")+
  theme_bw()


##1000 m2 richness
thousand_m$Year = as.factor(thousand_m$Year)
gamma <- thousand_m%>%
  dplyr::select(Year, Transect,  species)%>%
  #revalue(as.factor(Year), c("twentythree"="2023", "twentyfour"="2024"))
  mutate(Year = as.factor(Year))%>%
  pivot_wider(names_from = Year, values_from = species)%>%
  left_join(transect.info)#%>%
  #dplyr::rename(twentythree = 2023, )
  #dpylr::summarize(diff = 2024-)

gamma$diff <- gamma$'2024' - gamma$'2023'


ggplot(gamma, aes(Treatment, diff))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  geom_hline(yintercept = 0)+
  ylab("Change in gamma diversity (1 km2)")+
  theme_bw()



###distrbance covariates

met_disturb_2024 <- left_join(metrics, disturbance, by = c("Transect", "SoilVeg", "Treatment"))%>%
  subset(Year == "2024")


met_disturb_2024%>%
  subset(Treatment != "Reference")%>%
  ggplot(aes(perc_disturbance_2024, intercept, color = SoilVeg))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()


met_disturb_2024%>%
  subset(Treatment != "Reference")%>%
  ggplot(aes(perc_disturbance_2024, slope, color = SoilVeg))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()





