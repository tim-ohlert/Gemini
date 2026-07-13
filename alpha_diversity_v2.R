###Generic alpha diversity measures and indicator species analysis
###Modified to produce analyses for 2024, 2023, and change from 2023-2024

library(tidyverse)
library(plyr)
library(nlme)
library(emmeans)
library(vegan)
library(codyn)
library(indicspecies)
library(ggthemes)

modwhit.one <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/modwhit_clean_2024.csv")%>%
  subset(Quad_sz_m2 == 1)%>%
  subset(Cover_perc != "NA")

Mod.whit.spp <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/Mod-whit-spp.csv")
transect.info <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/transect-id-2024.csv")
disturbance <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/modwhit-disturbance.csv")
disturbance_quad <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/modwhit_disturbance.csv")
disturbance_quad$Quad_num <- as.character(disturbance_quad$Quad_num)

# #################################
# EVENNESS AND RICHNESS ANALYSIS
##################################

even_rich_all <- modwhit.one%>%
  unite("rep", c("Transect", "Quad_num"))%>%
  community_structure(time.var = "Year", abundance.var = "Cover_perc", 
                      replicate.var = "rep", metric = "Evar")%>%
  separate("rep", into = c("Transect", "Quad_num"))%>%
  left_join(disturbance_quad, by = c("Transect", "Quad_num"))%>%
  left_join(transect.info, by = "Transect")

even_rich_all$Treatment <- revalue(even_rich_all$Treatment, 
                                   c("Drive and Crush" = "Impact", "Reference" = "Control"))

#2024 Evenness
even_rich_2024 <- even_rich_all %>% subset(Year == 2024)

even_rich_2024 %>%
  ggplot(aes(Treatment, Evar, fill = Treatment))+
  facet_wrap(~SoilVeg)+
  scale_fill_manual(values = c( "#FFE793", "#4F93A7"))+
  geom_boxplot()+
  ylab("Evenness")+
  ylim(0,1)+
  xlab("")+
  theme_base()

ggsave("C:/Users/ohler/Dropbox/grants/Gemini/figures/evenness_2024.pdf",
       plot = last_plot(),
       device = "pdf",
       path = NULL,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in"),
       dpi = 600,
       limitsize = TRUE)

mod <- lme(Evar~Treatment*SoilVeg, random = ~1|Transect, data = subset(even_rich_2024, Evar != "NA"))
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

# 2024 Richness
even_rich_2024 %>%
  ggplot(aes(Treatment, richness, fill = Treatment))+
  facet_wrap(~SoilVeg)+
  scale_fill_manual(values = c( "#FFE793", "#4F93A7"))+
  geom_boxplot()+
  ylab("Richness")+
  xlab("")+
  ylim(0,20)+
  theme_base()

ggsave("C:/Users/ohler/Dropbox/grants/Gemini/figures/richness_2024.pdf",
       plot = last_plot(),
       device = "pdf",
       path = NULL,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in"),
       dpi = 600,
       limitsize = TRUE)

mod <- lme(richness~Treatment*SoilVeg, random = ~1|Transect, data = subset(even_rich_2024, Evar != "NA"))
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

#2023 Evenness
even_rich_2023 <- even_rich_all %>% subset(Year == 2023)

even_rich_2023 %>%
  ggplot(aes(Treatment, Evar, fill = Treatment))+
  facet_wrap(~SoilVeg)+
  scale_fill_manual(values = c( "#FFE793", "#4F93A7"))+
  geom_boxplot()+
  ylab("Evenness")+
  ylim(0,1)+
  xlab("")+
  theme_base()

ggsave("C:/Users/ohler/Dropbox/grants/Gemini/figures/evenness_2023.pdf",
       plot = last_plot(),
       device = "pdf",
       path = NULL,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in"),
       dpi = 600,
       limitsize = TRUE)

mod <- lme(Evar~Treatment*SoilVeg, random = ~1|Transect, data = subset(even_rich_2023, Evar != "NA"))
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

# 2023 Richness
even_rich_2023 %>%
  ggplot(aes(Treatment, richness, fill = Treatment))+
  facet_wrap(~SoilVeg)+
  scale_fill_manual(values = c( "#FFE793", "#4F93A7"))+
  geom_boxplot()+
  ylab("Richness")+
  xlab("")+
  ylim(0,20)+
  theme_base()

ggsave("C:/Users/ohler/Dropbox/grants/Gemini/figures/richness_2023.pdf",
       plot = last_plot(),
       device = "pdf",
       path = NULL,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in"),
       dpi = 600,
       limitsize = TRUE)

mod <- lme(richness~Treatment*SoilVeg, random = ~1|Transect, data = subset(even_rich_2023, Evar != "NA"))
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

#Richness and Evenness Change
even_change <- even_rich_all %>%
  dplyr::select(Year, Transect, Quad_num, Treatment, SoilVeg, Evar) %>%
  pivot_wider(names_from = Year, values_from = Evar) %>%
  mutate(Evar_change = `2024` - `2023`)

rich_change <- even_rich_all %>%
  dplyr::select(Year, Transect, Quad_num, Treatment, SoilVeg, richness) %>%
  pivot_wider(names_from = Year, values_from = richness) %>%
  mutate(richness_change = `2024` - `2023`)

even_change %>%
  ggplot(aes(Treatment, Evar_change, fill = Treatment))+
  facet_wrap(~SoilVeg)+
  scale_fill_manual(values = c( "#FFE793", "#4F93A7"))+
  geom_boxplot()+
  geom_hline(yintercept = 0, linetype = "dashed")+
  ylab("Change in Evenness")+
  xlab("")+
  theme_base()

ggsave("C:/Users/ohler/Dropbox/grants/Gemini/figures/evenness_change.pdf",
       plot = last_plot(),
       device = "pdf",
       path = NULL,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in"),
       dpi = 600,
       limitsize = TRUE)

mod <- lme(Evar_change~Treatment*SoilVeg, random = ~1|Transect, 
           data = even_change)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

rich_change %>%
  ggplot(aes(Treatment, richness_change, fill = Treatment))+
  facet_wrap(~SoilVeg)+
  scale_fill_manual(values = c( "#FFE793", "#4F93A7"))+
  geom_boxplot()+
  geom_hline(yintercept = 0, linetype = "dashed")+
  ylab("Change in Richness")+
  xlab("")+
  theme_base()

ggsave("C:/Users/ohler/Dropbox/grants/Gemini/figures/richness_change.pdf",
       plot = last_plot(),
       device = "pdf",
       path = NULL,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in"),
       dpi = 600,
       limitsize = TRUE)

mod <- lme(richness_change~Treatment*SoilVeg, random = ~1|Transect, 
           data = na.omit(rich_change))
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))



# ################################
# INDICATOR SPECIES ANALYSIS
#################################

wide_all <- modwhit.one%>%
  left_join(transect.info, by = "Transect")%>%
  dplyr::select(Year, Transect, Quad_num, Treatment, SoilVeg, Spp_code, Cover_perc)%>%
  pivot_wider( names_from = "Spp_code", values_from = "Cover_perc", values_fill = 0)

#2024 Indicator Species Analysis

#DEEPCREOSOTE 2024
x <- subset(wide_all, SoilVeg == "DeepCreosote" & Year == 2024)[,6:70]

indval <- multipatt(x
                    , subset(wide_all, SoilVeg == "DeepCreosote" & Year == 2024)$Treatment, 
                    control = how(nperm=999)) 
summary(indval)


#SHALLOWCREOSOTE 2024
x <- subset(wide_all, SoilVeg == "ShallowCreosote" & Year == 2024)[,6:70]

indval <- multipatt(x
                    , subset(wide_all, SoilVeg == "ShallowCreosote" & Year == 2024)$Treatment, 
                    control = how(nperm=999)) 
summary(indval)



#SILTYATRIPLEX 2024
x <- subset(wide_all, SoilVeg == "SiltyAtriplex" & Year == 2024)[,6:70]

indval <- multipatt(x
                    , subset(wide_all, SoilVeg == "SiltyAtriplex" & Year == 2024)$Treatment, 
                    control = how(nperm=999)) 
summary(indval)


# --- 2023 Indicator Species Analysis ---

#DEEPCREOSOTE 2023
x <- subset(wide_all, SoilVeg == "DeepCreosote" & Year == 2023)[,6:70]

indval <- multipatt(x
                    , subset(wide_all, SoilVeg == "DeepCreosote" & Year == 2023)$Treatment, 
                    control = how(nperm=999)) 
summary(indval)


#SHALLOWCREOSOTE 2023
x <- subset(wide_all, SoilVeg == "ShallowCreosote" & Year == 2023)[,6:70]

indval <- multipatt(x
                    , subset(wide_all, SoilVeg == "ShallowCreosote" & Year == 2023)$Treatment, 
                    control = how(nperm=999)) 
summary(indval)



#SILTYATRIPLEX 2023
x <- subset(wide_all, SoilVeg == "SiltyAtriplex" & Year == 2023)[,6:70]

indval <- multipatt(x
                    , subset(wide_all, SoilVeg == "SiltyAtriplex" & Year == 2023)$Treatment, 
                    control = how(nperm=999)) 
summary(indval)

