###Gneric alpha diversity measures and indicator species analysis

library(tidyverse)
library(plyr)
library(nlme)
library(emmeans)
library(vegan)
library(codyn)
library(indicspecies)





modwhit.one <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/modwhit_clean_2024.csv")%>%
              subset(Quad_sz_m2 == 1)%>%
              subset(Cover_perc != "NA")
Mod.whit.spp <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/Mod-whit-spp.csv")
transect.info <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/transect-id-2024.csv")
disturbance <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/modwhit-disturbance.csv")




shannon <- modwhit.one%>%
    unite("rep", c("Transect", "Quad_num" ))%>%
  community_diversity(time.var = "Year", abundance.var = "Cover_perc", replicate.var = "rep", metric = "Shannon")%>%
  separate("rep", into = c("Transect", "Quad_num"))%>%
  left_join(disturbance, "Transect")

shannon$Treatment <- revalue(shannon$Treatment, c("Drive and Crush" = "Impact", "Reference" = "Control" ))

shannon%>%
  subset(Year == 2024)%>%
ggplot( aes(Treatment, Shannon, color = Treatment))+
  facet_wrap(~SoilVeg)+
  scale_color_manual(values = c("black", "blue"))+
  geom_boxplot()+
  xlab("")+
  theme_bw()

mod <- lme(Shannon~Treatment*SoilVeg, random = ~1|Transect, data = subset(shannon, Year == 2024))
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))


even <- modwhit.one%>%
  unite("rep", c("Transect", "Quad_num" ))%>%
  community_structure(time.var = "Year", abundance.var = "Cover_perc", replicate.var = "rep", metric = "Evar")%>%
  separate("rep", into = c("Transect", "Quad_num"))%>%
  left_join(disturbance, "Transect")

even$Treatment <- revalue(even$Treatment, c("Drive and Crush" = "Impact", "Reference" = "Control" ))


even%>%
  subset(Year == 2024)%>%
ggplot( aes(Treatment, Evar, color = Treatment))+
  facet_wrap(~SoilVeg)+
  scale_color_manual(values = c("black", "blue"))+
  geom_boxplot()+
  xlab("")+
  theme_bw()

mod <- lme(Evar~Treatment*SoilVeg, random = ~1|Transect, data = subset(even, Year == 2024 & Evar != "NA"))
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))


even%>%
  subset(Year == 2024)%>%
  ggplot( aes(Treatment, richness, color = Treatment))+
  facet_wrap(~SoilVeg)+
  scale_color_manual(values = c("black", "blue"))+
  geom_boxplot()+
  xlab("")+
  ylim(0,20)+
  theme_bw()


mod <- lme(richness~Treatment*SoilVeg, random = ~1|Transect, data = subset(even, Year == 2024 & Evar != "NA"))
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))





###INDICATOR SPECIES ANALYSIS TRIAL


wide <- modwhit.one%>%
  left_join(transect.info, by = "Transect")%>%
  dplyr::select(Year, Transect, Quad_num, Treatment, SoilVeg, Spp_code, Cover_perc)%>%
  pivot_wider( names_from = "Spp_code", values_from = "Cover_perc", values_fill = 0)


#DEEPCREOSOTE
x <- subset(wide, SoilVeg == "DeepCreosote" & Year == 2024)[,6:70]

indval <- multipatt(x
                    , subset(wide, SoilVeg == "DeepCreosote" & Year == 2024)$Treatment, 
                    control = how(nperm=999)) 
summary(indval)


#SHALLOWCREOSOTE
x <- subset(wide, SoilVeg == "ShallowCreosote" & Year == 2024)[,6:70]

indval <- multipatt(x
                    , subset(wide, SoilVeg == "ShallowCreosote" & Year == 2024)$Treatment, 
                    control = how(nperm=999)) 
summary(indval)



#SILTYATRIPLEX
x <- subset(wide, SoilVeg == "SiltyAtriplex" & Year == 2024)[,6:70]

indval <- multipatt(x
                    , subset(wide, SoilVeg == "SiltyAtriplex" & Year == 2024)$Treatment, 
                    control = how(nperm=999)) 
summary(indval)




