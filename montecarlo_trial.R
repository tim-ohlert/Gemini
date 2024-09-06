
library(tidyverse)
library(plyr)
library(nlme)
library(emmeans)
library(vegan)
library(codyn)
library(tidyverse)
library(MASS) 

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

# Calculate correlation coefficient between mpg and wt

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
percentage_to_keep <- 1#0.95

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
  geom_histogram()

ggplot(mc, aes(Treatment, intercept, color = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
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
  geom_histogram()

ggplot(mc, aes(Treatment, slope, color = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  theme_bw()


mod <- lme(slope~Treatment, random = ~1|SoilVeg/Transect, data = mc)
summary(mod)

mod <- lme(slope~Treatment*SoilVeg, random = ~1|Transect, data = mc)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))
