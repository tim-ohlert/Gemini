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
visreg(mod)


#beta diversity (slope)
mod <- lme(slope~Treatment, random = ~1|SoilVeg, data = metrics)
summary(mod)

mod <- lm(slope~Treatment*SoilVeg, data = metrics)
summary(mod)
visreg(mod)


#gamma diversity (saturation point)


thousand_m <- rbind( thousand_m)%>%
  left_join(transect.info, by = "Transect")

mod <- lme(species~Treatment, random = ~1|SoilVeg, data = thousand_m)
summary(mod)

mod <- lm(species~Treatment*SoilVeg, data = thousand_m)
summary(mod)
visreg(mod)

##probably delete below
x <- plot_summ%>%
  subset(Year == "2024")%>%
  subset(Transect == "B14" & Year == "2024" & log_sp != "-Inf")

mod <- lm(species~log(Quad_sz_m2), data = x)
summary(mod)
newdata <- data.frame("Quad_sz_m2" = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100, 1000))#, 10000, 100000, 1000000, 10000000, 100000000,150000000,155000000, 1000000000,10000000000,100000000000,1000000000000,10000000000000))
y <- data.frame(predict(mod, newdata = newdata, interval = "prediction"))
plot(newdata$Quad_sz_m2, y$fit)

#%>%
  ggplot(aes(Quad_sz_m2, species, color = Treatment))+
  facet_wrap(~Transect)+
  geom_point()+
  geom_smooth(method = "glm", formula = y~log(x)


