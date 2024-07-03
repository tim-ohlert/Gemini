library(tidyverse)
library(plyr)


modwhit <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/2023_modwhit.csv")
transect.info <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/transect-id.csv")

all_but_thousand_m <- modwhit%>%
                  subset(Quad_sz_m2 != "1000")%>%
                  ddply(.(Transect, Quad_num, Quad_sz_m2), function(x)data.frame(
                    species = ifelse(x$Spp_code == "none", 0,
                                length(x$Spp_code)
                  )))

thousand_m <- modwhit%>%
           subset(Spp_code != "none") %>%
ddply(.(Transect), function(x)data.frame(
  species = length(unique(x$Spp_code))
))
thousand_m$Quad_num <- 14
thousand_m$Quad_sz_m2 <- 1000


plot_summ <- rbind(all_but_thousand_m, thousand_m)%>%
              unique()%>%
              left_join(transect.info, by = "Transect")

plot_summ$log_sp <- log(plot_summ$species)
plot_summ$log_area <- log(plot_summ$Quad_sz_m2)






