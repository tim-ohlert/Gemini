
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



#plot_summ%>%
#  subset(Year == "2024")%>%
#  ggplot(aes(log_area, log_sp, color = Treatment))+
#  facet_wrap(~SoilVeg)+
#  geom_point()+
#  geom_smooth(method = "lm")+
#  theme_bw()

#plot_summ%>%
#  subset(Year == "2024")%>%
#  ggplot(aes(Quad_sz_m2, species, color = Treatment))+
#  facet_wrap(~Transect)+
#  geom_point()+
#  geom_smooth(method = "glm", formula = y~log(x)#,
              #method.args = list(family = gaussian(link = 'log'))
#  )+
#  theme_bw()

#plot_summ%>%
#  subset(Year == "2024")%>%
#  ggplot(aes(log_area, log_sp, color = Treatment))+
#  facet_wrap(~Transect)+
#  geom_point()+
#  geom_smooth(method = "lm"#,
              #method.args = list(family = gaussian(link = 'log'))
#  )+
#  theme_bw()



#mod <- lme(log_sp~log_area*Treatment, random = ~1|SoilVeg/Transect, data = plot_summ%>%
#             subset(Year == "2024" & log_sp != "-Inf"))
#summary(mod)


#mod <- lm(species~log_area, data = plot_summ%>%
#            subset(Transect == "B14" & Year == "2024" & log_sp != "-Inf"))
#summary(mod)$coefficients[1]


#metrics <- plot_summ%>%
#  subset(Year == 2024 )%>%
#  ddply(.(Transect, SoilVeg, Treatment), function(x)data.frame(
#    intercept = summary(lm(log_sp~log_area, data = subset(x, log_sp != "-Inf")))$coefficients[1],
#    slope = summary(lm(log_sp~log_area, data = subset(x, log_sp != "-Inf")))$coefficients[2]
#  )) #add gamma diversity metric (saturation of curve)



#try monte carlo

plot_summ.1 <- plot_summ.1%>%
  unite(SoilVegTrt.col,SoilVeg:Treatment)

SoilVegTrt <- plot_summ.1%>%
  dplyr::select(SoilVegTrt.col)%>%
  unique()
master_trimmed_results <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("intercept", "slope", "SoilVegTrt"))

for (i in SoilVegTrt) {

# Calculate correlation coefficient between mpg and wt

plot_summ <- subset(plot_summ.1, log_sp != "-Inf" & SoilVegTrt.col == i)
correlation <- cor(plot_summ$log_sp, plot_summ$log_area)

sd_sp <- sd(plot_summ$log_sp, na.rm = TRUE)
sd_area <- sd(plot_summ$log_area, na.rm = TRUE)

#plot_summ%>%
#  subset(Year == "2024")%>%
#  group_by(SoilVeg, Treatment)%>%
#  dplyr::summarize(sd_sp = sd(log_sp, na.rm = TRUE), sd_area = sd(log_area, na.rm = TRUE))


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
num_replications <- 250

mc_results <- replicate(num_replications, {
  synthetic_data <- generate_synthetic_data()
  fit_model(synthetic_data)
}, simplify = FALSE)

# Convert results to a data frame
mc_results <- as.data.frame(do.call(rbind, mc_results))

# Print the results
print(mc_results)



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

# Plot the data with trimmed confidence intervals
#ggplot(plot_summ, aes(x = log_area, y = log_sp)) +
#  geom_point() +
#  geom_abline(data = trimmed_results, aes(intercept = `(Intercept)`, slope = log_area), color = "#f7aa58", alpha = 0.3) +
#  ggtitle("Monte Carlo Simulation of Linear Regression Confidence Intervals") +
#  labs(subtitle = "Gemini")


######intercept of trimmed results is alpha diversity intercept
#####column log_area is beta diversity slope

new_trimmed_results <- data.frame(intercept = trimmed_results$`(Intercept)`, slope = trimmed_results$log_area, SoilVegTrt.col = i)


master_trimmed_results <- rbind(master_trimmed_results, new_trimmed_results)

rm(new_trimmed_results)

}

