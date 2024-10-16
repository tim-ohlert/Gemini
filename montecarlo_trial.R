
library(tidyverse)
library(plyr)
library(nlme)
library(emmeans)
library(vegan)
library(codyn)
library(tidyverse)
library(MASS) 
library(boot)

modwhit <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/modwhit_clean_2024.csv")
Mod.whit.spp <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/Mod-whit-spp.csv")
transect.info <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/transect-id-2024.csv")
disturbance <- read.csv("C:/Users/ohler/Dropbox/grants/Gemini/modwhit-disturbance.csv")
disturbance$Treatment <- plyr::revalue(disturbance$Treatment, c("Reference" = "Control", "Drive and Crush" = "Impact"))
disturbance$SoilVeg <- plyr::revalue(disturbance$SoilVeg, c("DeepCreosote" = "Deep Creosote", "ShallowCreosote" = "Shallow Creosote", "SiltyAtriplex" = "Silty Saltbush"))

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
thousand_m$Quad_sz_m2 <- ifelse(thousand_m$Transect == "B14" & thousand_m$Year == 2024, 1250, 1000) #one large plot B14 got sampled at a larger area


four_thousand_m <- modwhit%>%
  left_join(transect.info, by = "Transect")%>%
  subset(Spp_code != "none") %>%
  ddply(.(Year, SoilVeg, Treatment), function(x)data.frame(
    species = length(unique(x$Spp_code))
  ))%>%
  subset(Year == 2024)

plot_summ.1 <- rbind(all_but_thousand_m, thousand_m)%>%
  unique()%>%
  left_join(transect.info, by = "Transect")%>%
  subset(Year == 2024)

plot_summ.1$log_sp <- log(plot_summ.1$species)
plot_summ.1$log_area <- log(plot_summ.1$Quad_sz_m2)

plot_summ.1 <- dplyr::mutate(plot_summ.1, Treatment = fct_recode(Treatment, "Impact" = "Drive and Crush", "Control" = "Reference"), SoilVeg = fct_recode(SoilVeg, "Silty Saltbush" = "SiltyAtriplex", "Shallow Creosote" = "ShallowCreosote", "Deep Creosote" = "DeepCreosote"))



#try monte carlo

plot_summ.1 <- plot_summ.1%>%
  unite(SoilVegTrt.col,c("Transect","SoilVeg","Treatment"))

SoilVegTrt <- plot_summ.1%>%
  dplyr::select(SoilVegTrt.col)
SoilVegTrt <- unique(SoilVegTrt$SoilVegTrt.col)

master_trimmed_results <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("intercept", "slope", "SoilVegTrt"))

for (i in 1:length(SoilVegTrt)) {

# Calculate correlation coefficient 

plot_summ <- subset(plot_summ.1, log_sp != "-Inf" & SoilVegTrt.col == SoilVegTrt[i])
correlation <- cor(plot_summ$log_sp, plot_summ$log_area)

sd_sp <- sd(plot_summ$log_sp, na.rm = TRUE)
sd_area <- sd(plot_summ$log_area, na.rm = TRUE)

# Define covariance matrix based on correlation coefficient
cov_matrix <- matrix(c(sd_sp^2, correlation * sd_area * sd_sp,
                       correlation * sd_sp * sd_area, sd_area^2), nrow = 2)

# Define function to generate synthetic data
generate_synthetic_data <- function() {
  synthetic_data <- mvrnorm(n = nrow(plot_summ), mu = c(mean(plot_summ$log_sp), mean(plot_summ$log_area)), Sigma = cov_matrix)
  synthetic_data <- as.data.frame(synthetic_data)
  colnames(synthetic_data) <- c("log_sp", "log_area")
  return(synthetic_data)
}

# Define function to fit linear regression model
fit_model <- function(data) {
  fit <- lm(log_sp ~ log_area, data = data)
  return(coef(fit))
}

# Perform Monte Carlo simulation
num_replications <- 1000

mc_results <- replicate(num_replications, {
  synthetic_data <- generate_synthetic_data()
  fit_model(synthetic_data)
}, simplify = FALSE)

# Convert results to a data frame
mc_results <- as.data.frame(do.call(rbind, mc_results))

# Print the results
#print(mc_results)

# Calculate the centroid of the cluster of regression lines
centroid <- colMeans(mc_results)

# Calculate the Euclidean distance of each line from the centroid
distances <- sqrt(rowSums((mc_results - centroid)^2))

# Specify the percentage of lines to keep (e.g., 95%)
percentage_to_keep <- 0.95

# Determine the number of lines to keep based on the specified percentage
num_lines_to_keep <- ceiling(length(distances) * percentage_to_keep)

# Find the indices of the lines closest to the centroid
indices_to_keep <- order(distances)[1:num_lines_to_keep]

# Filter out lines that are farthest from the centroid
trimmed_results <- mc_results[indices_to_keep, ]


######intercept of trimmed results is alpha diversity intercept
#####column log_area is beta diversity slope

new_trimmed_results <- data.frame(intercept = trimmed_results$`(Intercept)`, slope = trimmed_results$log_area, SoilVegTrt.col = SoilVegTrt[i])


master_trimmed_results <- rbind(master_trimmed_results, new_trimmed_results)

rm(new_trimmed_results)

}

mc <- separate(master_trimmed_results, SoilVegTrt.col, into = c("Transect", "SoilVeg", "Treatment"), sep = "_")




##intercept
ggplot(mc, aes(intercept))+
  facet_grid(Treatment~SoilVeg)+
  geom_histogram()+
  xlab("log(c)")+
  xlim(0,3)+
  theme_bw()

mean_vals <- mc%>%
  group_by(SoilVeg, Treatment)%>%
  dplyr::summarise(intercept = mean(intercept), slope = mean(slope))

ggplot(mc, aes())+
  facet_grid(Treatment~SoilVeg)+
  geom_abline(aes(intercept = intercept, slope = slope), alpha = 0.01)+
  geom_abline(data=mean_vals, aes(intercept = intercept, slope = slope), alpha = 1, color = "dodgerblue")+
  ylim(0,5)+
  xlim(0,15)+
  xlab("log(area)")+
  ylab("log(species)")+
  theme_bw()

ggplot(mc, aes(Treatment, intercept, color = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  scale_color_manual(values = c("black", "blue"))+
  xlab("")+
  ylab("log(c)")+
  ylim(0,3)+
  theme_bw()

mod <- lme(intercept~Treatment, random = ~1|SoilVeg/Transect, data = mc)
summary(mod)

mod <- lme(intercept~Treatment*SoilVeg, random = ~1|Transect, data = mc)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

  

##slope
ggplot(mc, aes(slope))+
  facet_grid(Treatment~SoilVeg)+
  geom_histogram()+
  xlab("z")+
  xlim(-0.35,0.75)+
  theme_bw()

ggplot(mc, aes(Treatment, slope, color = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  scale_color_manual(values = c("black", "blue"))+
  xlab("")+
  theme_bw()


mod <- lme(slope~Treatment, random = ~1|SoilVeg/Transect, data = mc)
summary(mod)

mod <- lme(slope~Treatment*SoilVeg, random = ~1|Transect, data = mc)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))




#################################
#####bootstrapping gamma diversity

four_thousand_m <- modwhit%>%
  left_join(transect.info, by = "Transect")%>%
  subset(Spp_code != "none") %>%
  subset(Year == 2024)%>%
  dplyr::select(Transect, Spp_code, SoilVeg, Treatment)%>%
  dplyr::select( Spp_code, SoilVeg, Treatment)%>%
  unite(SoilVegTrt.col,c("SoilVeg","Treatment"))
 # unique()
  #ddply(.(SoilVeg, Treatment), function(x)data.frame(
  #  species = length(unique(x$Spp_code))
  #))
  

SoilVegTrt <- four_thousand_m%>%
  dplyr::select(SoilVegTrt.col)
SoilVegTrt <- unique(SoilVegTrt$SoilVegTrt.col)

master_gamma_results <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("SoilVegTrt.col", "sr"))


for (i in 1:length(SoilVegTrt)) {

sr_boot.temp <-  replicate(1000,{
  d <- subset(four_thousand_m, SoilVegTrt.col == SoilVegTrt[i])
  
  temp <- sample(x = d$Spp_code, size = length(d$Spp_code), replace = TRUE)
    temp.sr <- as.numeric(length(unique(temp)))
  } )

new_results <- data.frame(SoilVegTrt.col = SoilVegTrt[i], sr = sr_boot.temp)

master_gamma_results <- rbind(master_gamma_results, new_results)

rm(new_results)
rm(sr_boot.temp)
}
  


master_gamma_results <- master_gamma_results%>%
  separate(SoilVegTrt.col, c("SoilVeg","Treatment"), sep = "_")

master_gamma_results$Treatment <- revalue(master_gamma_results$Treatment, c("Drive and Crush" = "Impact", "Reference" = "Control" ))

ggplot(master_gamma_results, aes(Treatment, sr, color = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  scale_color_manual(values = c("black", "blue"))+
  ylim(0,60)+
  xlab("")+
  ylab("Richness (4,000 m2)")+
  theme_bw()

mod <- lm(sr~Treatment*SoilVeg, data = master_gamma_results)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

  
######################
############
  

###distrbance covariates

met_disturb_2024 <- left_join(mc, disturbance, by = c("Transect", "SoilVeg", "Treatment"))
#revalue(met_disturb_2024$Treatment, c("Control" = "Reference", "Drive and Crush" = "Impact"))


met_disturb_2024%>%
  subset(Treatment == "Impact")%>%
  ggplot(aes(perc_disturbance_2024, intercept, color = SoilVeg))+
  #facet_wrap(~SoilVeg)+
  geom_point(alpha = 0.01)+
  geom_smooth(method = "lm")+
  ylab("Intercept (log(c)")+
  xlab("Percent disturbance")+
  theme_bw()


met_disturb_2024%>%
  subset(Treatment != "Control")%>%
  ggplot(aes(perc_disturbance_2024, slope, color = SoilVeg))+
  geom_point(alpha = 0.01)+
  geom_smooth(method = "lm")+
  ylab("Slope (z)")+
  xlab("Percent disturbance")+
  theme_bw()

