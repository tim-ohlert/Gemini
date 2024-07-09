library(tidyverse)
library(plyr)
library(nlme)

modwhit <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/modwhit_clean_2024.csv")
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
    
mod <- lme(log_sp~log_area*Treatment, random = ~1|SoilVeg/Transect, data = plot_summ%>%
             subset(Year == "2024" & log_sp != "-Inf"))
summary(mod)


mod <- lm(species~log_area, data = plot_summ%>%
             subset(Transect == "B14" & Year == "2024" & log_sp != "-Inf"))
summary(mod)$coefficients[1]

#alpha diversity (intercept)

metrics <- plot_summ%>%
  subset(Year == 2024 )%>%
  ddply(.(Transect, SoilVeg, Treatment), function(x)data.frame(
    intercept = summary(lm(log_sp~log_area, data = subset(x, log_sp != "-Inf")))$coefficients[1]#,
    slope = summary(lm(log_sp~log_area, data = subset(x, log_sp != "-Inf")))$coefficients[2]
  ))


