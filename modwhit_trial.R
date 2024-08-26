library(tidyverse)
library(plyr)
library(nlme)
library(emmeans)
library(vegan)
library(codyn)

modwhit <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/modwhit_clean_2024.csv")
Mod.whit.spp <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/Mod-whit-spp.csv")
transect.info <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/transect-id-2024.csv")

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


thousand_m <- rbind( thousand_m)%>%
  left_join(transect.info, by = "Transect")%>%
  dplyr::mutate( Treatment = fct_recode(Treatment, "Impact" = "Drive and Crush", "Control" = "Reference"), SoilVeg = fct_recode(SoilVeg, "Silty Saltbush" = "SiltyAtriplex", "Shallow Creosote" = "ShallowCreosote", "Deep Creosote" = "DeepCreosote"))


mod <- lme(species~Treatment, random = ~1|SoilVeg, data = thousand_m)
summary(mod)


mod <- lm(species~Treatment*SoilVeg, data = thousand_m)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))
#visreg(mod)


ggplot(thousand_m, aes(Treatment, species))+
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


thousand_m <- rbind( thousand_m)%>%
  left_join(transect.info, by = "Transect")%>%
  dplyr::mutate( Treatment = fct_recode(Treatment, "Impact" = "Drive and Crush", "Control" = "Reference"), SoilVeg = fct_recode(SoilVeg, "Silty Saltbush" = "SiltyAtriplex", "Shallow Creosote" = "ShallowCreosote", "Deep Creosote" = "DeepCreosote"))


mod <- lme(species~Treatment, random = ~1|SoilVeg, data = thousand_m)
summary(mod)


mod <- lm(species~Treatment*SoilVeg, data = thousand_m)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))
#visreg(mod)


ggplot(thousand_m, aes(Treatment, species))+
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


thousand_m <- rbind( thousand_m)%>%
  left_join(transect.info, by = "Transect")%>%
  dplyr::mutate( Treatment = fct_recode(Treatment, "Impact" = "Drive and Crush", "Control" = "Reference"), SoilVeg = fct_recode(SoilVeg, "Silty Saltbush" = "SiltyAtriplex", "Shallow Creosote" = "ShallowCreosote", "Deep Creosote" = "DeepCreosote"))


mod <- lme(species~Treatment, random = ~1|SoilVeg, data = thousand_m)
summary(mod)


mod <- lm(species~Treatment*SoilVeg, data = thousand_m)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))
#visreg(mod)


ggplot(thousand_m, aes(Treatment, species))+
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
  subset(Quad_sz_m2 == 1)%>%
  subset(Year == 2024)%>%
  dplyr::select(Transect, Quad_num,Spp_code, Cover_perc)%>%
  pivot_wider(names_from = "Spp_code", values_from = "Cover_perc", values_fill = 0)%>%
  pivot_longer(3:63, names_to = "species", values_to = "cover")%>%
  left_join(transect.info, by = "Transect")%>%
  group_by(SoilVeg, Treatment, species)%>%
  dplyr::summarize(cover = mean(cover))

rank_abun%>%
  subset(cover > 0.5)%>%
ggplot( aes(species, cover))+
  facet_grid(SoilVeg~Treatment)+
  geom_bar(stat = "identity")




