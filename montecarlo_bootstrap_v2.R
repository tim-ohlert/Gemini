library(tidyverse)
library(plyr)
library(nlme)
library(emmeans)
library(vegan)
library(codyn)
library(tidyverse)
library(MASS) 
library(boot)
library(ggthemes)
library(MuMIn)

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
thousand_m$Quad_sz_m2 <- ifelse(thousand_m$Transect == "B14" & thousand_m$Year == 2024, 1250, 1000)

# ###################################
# 2024 ALPHA/BETA DIVERSITY ANALYSIS
# ######################################
plot_summ.1 <- rbind(all_but_thousand_m, thousand_m)%>%
  unique()%>%
  left_join(transect.info, by = "Transect")%>%
  subset(Year == 2024)

plot_summ.1$log_sp <- log(plot_summ.1$species)
plot_summ.1$log_area <- log(plot_summ.1$Quad_sz_m2)

plot_summ.1 <- dplyr::mutate(plot_summ.1, Treatment = fct_recode(Treatment, "Impact" = "Drive and Crush", "Control" = "Reference"), SoilVeg = fct_recode(SoilVeg, "Silty Saltbush" = "SiltyAtriplex", "Shallow Creosote" = "ShallowCreosote", "Deep Creosote" = "DeepCreosote"))

plot_summ.1 <- plot_summ.1%>%
  unite(SoilVegTrt.col,c("Transect","SoilVeg","Treatment"))

SoilVegTrt <- plot_summ.1%>%
  dplyr::select(SoilVegTrt.col)
SoilVegTrt <- unique(SoilVegTrt$SoilVegTrt.col)

master_trimmed_results_2024 <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("intercept", "slope", "SoilVegTrt"))

for (i in 1:length(SoilVegTrt)) {
  
  plot_summ <- subset(plot_summ.1, log_sp != "-Inf" & SoilVegTrt.col == SoilVegTrt[i])
  correlation <- cor(plot_summ$log_sp, plot_summ$log_area)
  
  sd_sp <- sd(plot_summ$log_sp, na.rm = TRUE)
  sd_area <- sd(plot_summ$log_area, na.rm = TRUE)
  
  cov_matrix <- matrix(c(sd_sp^2, correlation * sd_area * sd_sp,
                         correlation * sd_sp * sd_area, sd_area^2), nrow = 2)
  
  generate_synthetic_data <- function() {
    synthetic_data <- mvrnorm(n = nrow(plot_summ), mu = c(mean(plot_summ$log_sp), mean(plot_summ$log_area)), Sigma = cov_matrix)
    synthetic_data <- as.data.frame(synthetic_data)
    colnames(synthetic_data) <- c("log_sp", "log_area")
    return(synthetic_data)
  }
  
  fit_model <- function(data) {
    fit <- lm(log_sp ~ log_area, data = data)
    return(coef(fit))
  }
  
  num_replications <- 1000
  
  mc_results <- replicate(num_replications, {
    synthetic_data <- generate_synthetic_data()
    fit_model(synthetic_data)
  }, simplify = FALSE)
  
  mc_results <- as.data.frame(do.call(rbind, mc_results))
  
  centroid <- colMeans(mc_results)
  distances <- sqrt(rowSums((mc_results - centroid)^2))
  
  percentage_to_keep <- 0.95
  num_lines_to_keep <- ceiling(length(distances) * percentage_to_keep)
  indices_to_keep <- order(distances)[1:num_lines_to_keep]
  trimmed_results <- mc_results[indices_to_keep, ]
  
  new_trimmed_results <- data.frame(intercept = trimmed_results$`(Intercept)`, slope = trimmed_results$log_area, SoilVegTrt.col = SoilVegTrt[i])
  
  master_trimmed_results_2024 <- rbind(master_trimmed_results_2024, new_trimmed_results)
  
  rm(new_trimmed_results)
}

mc_2024 <- separate(master_trimmed_results_2024, SoilVegTrt.col, into = c("Transect", "SoilVeg", "Treatment"), sep = "_")

ggplot(mc_2024, aes(intercept))+
  facet_grid(Treatment~SoilVeg)+
  geom_histogram()+
  xlab("log(c)")+
  xlim(0,3)+
  theme_base()

mean_vals <- mc_2024%>%
  group_by(SoilVeg, Treatment)%>%
  dplyr::summarise(intercept = mean(intercept), slope = mean(slope))

ggplot(mc_2024, aes())+
  facet_grid(Treatment~SoilVeg)+
  geom_abline(aes(intercept = intercept, slope = slope), alpha = 0.01)+
  geom_abline(data=mean_vals, aes(intercept = intercept, slope = slope), alpha = 1, color = "dodgerblue")+
  ylim(0,5)+
  xlim(0,15)+
  xlab("log(area)")+
  ylab("log(species)")+
  theme_base()

ggplot(mc_2024, aes(Treatment, intercept, color = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  scale_color_manual(values = c("black", "blue"))+
  xlab("")+
  ylab("Alpha diversity [log(c)]")+
  ylim(0,3)+
  theme_base()

mod <- lme(intercept~Treatment*SoilVeg, random = ~1|Transect, data = mc_2024)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

ggplot(mc_2024, aes(Treatment, slope, fill = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  scale_fill_manual(values = c( "#FFE793", "#4F93A7"))+
  xlab("")+
  ylab("Beta diversity (z)")+
  theme_base()

ggsave("C:/Users/ohler/Dropbox/grants/Gemini/figures/betadiv_2024.pdf",
       plot = last_plot(),
       device = "pdf",
       path = NULL,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in"),
       dpi = 600,
       limitsize = TRUE)

mod <- lme(slope~Treatment*SoilVeg, random = ~1|Transect, data = mc_2024)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

# ============================================================================
# 2023 ALPHA/BETA DIVERSITY ANALYSIS
# ============================================================================

plot_summ.1_2023 <- rbind(all_but_thousand_m, thousand_m)%>%
  unique()%>%
  left_join(transect.info, by = "Transect")%>%
  subset(Year == 2023)

plot_summ.1_2023$log_sp <- log(plot_summ.1_2023$species)
plot_summ.1_2023$log_area <- log(plot_summ.1_2023$Quad_sz_m2)

plot_summ.1_2023 <- dplyr::mutate(plot_summ.1_2023, Treatment = fct_recode(Treatment, "Impact" = "Drive and Crush", "Control" = "Reference"), SoilVeg = fct_recode(SoilVeg, "Silty Saltbush" = "SiltyAtriplex", "Shallow Creosote" = "ShallowCreosote", "Deep Creosote" = "DeepCreosote"))

plot_summ.1_2023 <- plot_summ.1_2023%>%
  unite(SoilVegTrt.col,c("Transect","SoilVeg","Treatment"))

SoilVegTrt_2023 <- plot_summ.1_2023%>%
  dplyr::select(SoilVegTrt.col)
SoilVegTrt_2023 <- unique(SoilVegTrt_2023$SoilVegTrt.col)

master_trimmed_results_2023 <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("intercept", "slope", "SoilVegTrt"))

for (i in 1:length(SoilVegTrt_2023)) {
  
  plot_summ <- subset(plot_summ.1_2023, log_sp != "-Inf" & SoilVegTrt.col == SoilVegTrt_2023[i])
  correlation <- cor(plot_summ$log_sp, plot_summ$log_area)
  
  sd_sp <- sd(plot_summ$log_sp, na.rm = TRUE)
  sd_area <- sd(plot_summ$log_area, na.rm = TRUE)
  
  cov_matrix <- matrix(c(sd_sp^2, correlation * sd_area * sd_sp,
                         correlation * sd_sp * sd_area, sd_area^2), nrow = 2)
  
  generate_synthetic_data <- function() {
    synthetic_data <- mvrnorm(n = nrow(plot_summ), mu = c(mean(plot_summ$log_sp), mean(plot_summ$log_area)), Sigma = cov_matrix)
    synthetic_data <- as.data.frame(synthetic_data)
    colnames(synthetic_data) <- c("log_sp", "log_area")
    return(synthetic_data)
  }
  
  fit_model <- function(data) {
    fit <- lm(log_sp ~ log_area, data = data)
    return(coef(fit))
  }
  
  num_replications <- 1000
  
  mc_results <- replicate(num_replications, {
    synthetic_data <- generate_synthetic_data()
    fit_model(synthetic_data)
  }, simplify = FALSE)
  
  mc_results <- as.data.frame(do.call(rbind, mc_results))
  
  centroid <- colMeans(mc_results)
  distances <- sqrt(rowSums((mc_results - centroid)^2))
  
  percentage_to_keep <- 0.95
  num_lines_to_keep <- ceiling(length(distances) * percentage_to_keep)
  indices_to_keep <- order(distances)[1:num_lines_to_keep]
  trimmed_results <- mc_results[indices_to_keep, ]
  
  new_trimmed_results <- data.frame(intercept = trimmed_results$`(Intercept)`, slope = trimmed_results$log_area, SoilVegTrt.col = SoilVegTrt_2023[i])
  
  master_trimmed_results_2023 <- rbind(master_trimmed_results_2023, new_trimmed_results)
  
  rm(new_trimmed_results)
}

mc_2023 <- separate(master_trimmed_results_2023, SoilVegTrt.col, into = c("Transect", "SoilVeg", "Treatment"), sep = "_")

ggplot(mc_2023, aes(intercept))+
  facet_grid(Treatment~SoilVeg)+
  geom_histogram()+
  xlab("log(c)")+
  xlim(0,3)+
  theme_base()

mean_vals_2023 <- mc_2023%>%
  group_by(SoilVeg, Treatment)%>%
  dplyr::summarise(intercept = mean(intercept), slope = mean(slope))

ggplot(mc_2023, aes())+
  facet_grid(Treatment~SoilVeg)+
  geom_abline(aes(intercept = intercept, slope = slope), alpha = 0.01)+
  geom_abline(data=mean_vals_2023, aes(intercept = intercept, slope = slope), alpha = 1, color = "dodgerblue")+
  ylim(0,5)+
  xlim(0,15)+
  xlab("log(area)")+
  ylab("log(species)")+
  theme_base()

ggplot(mc_2023, aes(Treatment, intercept, color = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  scale_color_manual(values = c("black", "blue"))+
  xlab("")+
  ylab("Alpha diversity [log(c)]")+
  ylim(0,3)+
  theme_base()

mod <- lme(intercept~Treatment*SoilVeg, random = ~1|Transect, data = mc_2023)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

ggplot(mc_2023, aes(Treatment, slope, fill = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  scale_fill_manual(values = c( "#FFE793", "#4F93A7"))+
  xlab("")+
  ylab("Beta diversity (z)")+
  theme_base()

ggsave("C:/Users/ohler/Dropbox/grants/Gemini/figures/betadiv_2023.pdf",
       plot = last_plot(),
       device = "pdf",
       path = NULL,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in"),
       dpi = 600,
       limitsize = TRUE)

mod <- lme(slope~Treatment*SoilVeg, random = ~1|Transect, data = mc_2023)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

# ============================================================================
# ALPHA/BETA DIVERSITY CHANGE (2024-2023)
# ============================================================================

mc_change <- mc_2024 %>%
  dplyr::group_by(Transect, SoilVeg, Treatment) %>%
  dplyr::mutate(replicate = row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::rename(intercept_2024 = intercept, slope_2024 = slope) %>%
  left_join(
    mc_2023 %>%
      dplyr::group_by(Transect, SoilVeg, Treatment) %>%
      dplyr::mutate(replicate = row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::select(Transect, SoilVeg, Treatment, replicate, 
                    intercept_2023 = intercept, slope_2023 = slope),
    by = c("Transect", "SoilVeg", "Treatment", "replicate")
  ) %>%
  dplyr::mutate(intercept_change = intercept_2024 - intercept_2023,
                slope_change = slope_2024 - slope_2023) %>%
  na.omit()

ggplot(mc_change, aes(Treatment, intercept_change, color = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  scale_color_manual(values = c("black", "blue"))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  xlab("")+
  ylab("Change in alpha diversity [log(c)]")+
  theme_base()

mod <- lme(intercept_change~Treatment*SoilVeg, random = ~1|Transect, data = mc_change)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

ggplot(mc_change, aes(Treatment, slope_change, fill = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  scale_fill_manual(values = c( "#FFE793", "#4F93A7"))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  xlab("")+
  ylab("Change in beta diversity (z)")+
  theme_base()

ggsave("C:/Users/ohler/Dropbox/grants/Gemini/figures/betadiv_change.pdf",
       plot = last_plot(),
       device = "pdf",
       path = NULL,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in"),
       dpi = 600,
       limitsize = TRUE)

mod <- lme(slope_change~Treatment*SoilVeg, random = ~1|Transect, data = mc_change)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

# ============================================================================
# 2024 GAMMA DIVERSITY BOOTSTRAPPING (4000m)
# ============================================================================

four_thousand_m_2024 <- modwhit%>%
  left_join(transect.info, by = "Transect")%>%
  subset(Spp_code != "none") %>%
  subset(Year == 2024)%>%
  dplyr::select(Transect, Spp_code, SoilVeg, Treatment)%>%
  unique()%>%
  dplyr::select(Spp_code, SoilVeg, Treatment)%>%
  unite(SoilVegTrt.col,c("SoilVeg","Treatment"))

SoilVegTrt_4000_2024 <- four_thousand_m_2024%>%
  dplyr::select(SoilVegTrt.col)
SoilVegTrt_4000_2024 <- unique(SoilVegTrt_4000_2024$SoilVegTrt.col)

master_gamma_results_2024_4000 <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("SoilVegTrt.col", "sr"))

for (i in 1:length(SoilVegTrt_4000_2024)) {
  
  sr_boot.temp <-  replicate(1000,{
    d <- subset(four_thousand_m_2024, SoilVegTrt.col == SoilVegTrt_4000_2024[i])
    temp <- sample(x = d$Spp_code, size = length(d$Spp_code), replace = TRUE)
    temp.sr <- as.numeric(length(unique(temp)))
  } )
  
  new_results <- data.frame(SoilVegTrt.col = SoilVegTrt_4000_2024[i], sr = sr_boot.temp)
  
  master_gamma_results_2024_4000 <- rbind(master_gamma_results_2024_4000, new_results)
  
  rm(new_results)
  rm(sr_boot.temp)
}

master_gamma_results_2024_4000 <- master_gamma_results_2024_4000%>%
  separate(SoilVegTrt.col, c("SoilVeg","Treatment"), sep = "_")

master_gamma_results_2024_4000$Treatment <- revalue(master_gamma_results_2024_4000$Treatment, c("Drive and Crush" = "Impact", "Reference" = "Control"))

ggplot(master_gamma_results_2024_4000, aes(Treatment, sr, fill = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  scale_fill_manual(values = c( "#FFE793", "#4F93A7"))+
  ylim(0,60)+
  xlab("")+
  ylab("Gamma diversity (4,000 m2 richness)")+
  theme_base()

ggsave("C:/Users/ohler/Dropbox/grants/Gemini/figures/gammadiv_4000_2024.pdf",
       plot = last_plot(),
       device = "pdf",
       path = NULL,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in"),
       dpi = 600,
       limitsize = TRUE)

mod <- lm(sr~Treatment*SoilVeg, data = master_gamma_results_2024_4000)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

# ============================================================================
# 2023 GAMMA DIVERSITY BOOTSTRAPPING (4000m)
# ============================================================================

four_thousand_m_2023 <- modwhit%>%
  left_join(transect.info, by = "Transect")%>%
  subset(Spp_code != "none") %>%
  subset(Year == 2023)%>%
  dplyr::select(Transect, Spp_code, SoilVeg, Treatment)%>%
  unique()%>%
  dplyr::select(Spp_code, SoilVeg, Treatment)%>%
  unite(SoilVegTrt.col,c("SoilVeg","Treatment"))

SoilVegTrt_4000_2023 <- four_thousand_m_2023%>%
  dplyr::select(SoilVegTrt.col)
SoilVegTrt_4000_2023 <- unique(SoilVegTrt_4000_2023$SoilVegTrt.col)

master_gamma_results_2023_4000 <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("SoilVegTrt.col", "sr"))

for (i in 1:length(SoilVegTrt_4000_2023)) {
  
  sr_boot.temp <-  replicate(1000,{
    d <- subset(four_thousand_m_2023, SoilVegTrt.col == SoilVegTrt_4000_2023[i])
    temp <- sample(x = d$Spp_code, size = length(d$Spp_code), replace = TRUE)
    temp.sr <- as.numeric(length(unique(temp)))
  } )
  
  new_results <- data.frame(SoilVegTrt.col = SoilVegTrt_4000_2023[i], sr = sr_boot.temp)
  
  master_gamma_results_2023_4000 <- rbind(master_gamma_results_2023_4000, new_results)
  
  rm(new_results)
  rm(sr_boot.temp)
}

master_gamma_results_2023_4000 <- master_gamma_results_2023_4000%>%
  separate(SoilVegTrt.col, c("SoilVeg","Treatment"), sep = "_")

master_gamma_results_2023_4000$Treatment <- revalue(master_gamma_results_2023_4000$Treatment, c("Drive and Crush" = "Impact", "Reference" = "Control"))

ggplot(master_gamma_results_2023_4000, aes(Treatment, sr, fill = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  scale_fill_manual(values = c( "#FFE793", "#4F93A7"))+
  ylim(0,60)+
  xlab("")+
  ylab("Gamma diversity (4,000 m2 richness)")+
  theme_base()

ggsave("C:/Users/ohler/Dropbox/grants/Gemini/figures/gammadiv_4000_2023.pdf",
       plot = last_plot(),
       device = "pdf",
       path = NULL,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in"),
       dpi = 600,
       limitsize = TRUE)

mod <- lm(sr~Treatment*SoilVeg, data = master_gamma_results_2023_4000)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

# ============================================================================
# 2024 GAMMA DIVERSITY BOOTSTRAPPING (1000m)
# ============================================================================

one_thousand_m_2024 <- modwhit%>%
  left_join(transect.info, by = "Transect")%>%
  subset(Spp_code != "none") %>%
  subset(Year == 2024)%>%
  dplyr::select(Transect, Spp_code, SoilVeg, Treatment)%>%
  unite(SoilVegTrtTransect.col,c("SoilVeg","Treatment", "Transect"))

SoilVegTrtTransect_2024 <- one_thousand_m_2024%>%
  dplyr::select(SoilVegTrtTransect.col)
SoilVegTrtTransect_2024 <- unique(SoilVegTrtTransect_2024$SoilVegTrtTransect.col)

master_gamma_results_2024 <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("SoilVegTrtTransect.col", "sr"))

for (i in 1:length(SoilVegTrtTransect_2024)) {
  
  sr_boot.temp <-  replicate(1000,{
    d <- subset(one_thousand_m_2024, SoilVegTrtTransect.col == SoilVegTrtTransect_2024[i])
    temp <- sample(x = d$Spp_code, size = length(d$Spp_code), replace = TRUE)
    temp.sr <- as.numeric(length(unique(temp)))
  } )
  
  new_results <- data.frame(SoilVegTrtTransect.col = SoilVegTrtTransect_2024[i], sr = sr_boot.temp)
  
  master_gamma_results_2024 <- rbind(master_gamma_results_2024, new_results)
  
  rm(new_results)
  rm(sr_boot.temp)
}

master_gamma_results_2024 <- master_gamma_results_2024%>%
  separate(SoilVegTrtTransect.col, c("SoilVeg","Treatment", "Transect"), sep = "_")

master_gamma_results_2024$Treatment <- revalue(master_gamma_results_2024$Treatment, c("Drive and Crush" = "Impact", "Reference" = "Control"))

ggplot(master_gamma_results_2024, aes(Treatment, sr, fill = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  scale_fill_manual(values = c( "#FFE793", "#4F93A7"))+
  ylim(0,45)+
  xlab("")+
  ylab("Gamma diversity (1,000 m2 richness)")+
  theme_base()

ggsave("C:/Users/ohler/Dropbox/grants/Gemini/figures/gammadiv_2024.pdf",
       plot = last_plot(),
       device = "pdf",
       path = NULL,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in"),
       dpi = 600,
       limitsize = TRUE)

mod <- lme(sr~Treatment*SoilVeg, random = ~1|Transect, data = master_gamma_results_2024)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

# ============================================================================
# 2023 GAMMA DIVERSITY BOOTSTRAPPING (1000m)
# ============================================================================

one_thousand_m_2023 <- modwhit%>%
  left_join(transect.info, by = "Transect")%>%
  subset(Spp_code != "none") %>%
  subset(Year == 2023)%>%
  dplyr::select(Transect, Spp_code, SoilVeg, Treatment)%>%
  unite(SoilVegTrtTransect.col,c("SoilVeg","Treatment", "Transect"))

SoilVegTrtTransect_2023 <- one_thousand_m_2023%>%
  dplyr::select(SoilVegTrtTransect.col)
SoilVegTrtTransect_2023 <- unique(SoilVegTrtTransect_2023$SoilVegTrtTransect.col)

master_gamma_results_2023 <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("SoilVegTrtTransect.col", "sr"))

for (i in 1:length(SoilVegTrtTransect_2023)) {
  
  sr_boot.temp <-  replicate(1000,{
    d <- subset(one_thousand_m_2023, SoilVegTrtTransect.col == SoilVegTrtTransect_2023[i])
    temp <- sample(x = d$Spp_code, size = length(d$Spp_code), replace = TRUE)
    temp.sr <- as.numeric(length(unique(temp)))
  } )
  
  new_results <- data.frame(SoilVegTrtTransect.col = SoilVegTrtTransect_2023[i], sr = sr_boot.temp)
  
  master_gamma_results_2023 <- rbind(master_gamma_results_2023, new_results)
  
  rm(new_results)
  rm(sr_boot.temp)
}

master_gamma_results_2023 <- master_gamma_results_2023%>%
  separate(SoilVegTrtTransect.col, c("SoilVeg","Treatment", "Transect"), sep = "_")

master_gamma_results_2023$Treatment <- revalue(master_gamma_results_2023$Treatment, c("Drive and Crush" = "Impact", "Reference" = "Control"))

ggplot(master_gamma_results_2023, aes(Treatment, sr, fill = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  scale_fill_manual(values = c( "#FFE793", "#4F93A7"))+
  ylim(0,45)+
  xlab("")+
  ylab("Gamma diversity (1,000 m2 richness)")+
  theme_base()

ggsave("C:/Users/ohler/Dropbox/grants/Gemini/figures/gammadiv_2023.pdf",
       plot = last_plot(),
       device = "pdf",
       path = NULL,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in"),
       dpi = 600,
       limitsize = TRUE)

mod <- lme(sr~Treatment*SoilVeg, random = ~1|Transect, data = master_gamma_results_2023)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

# ============================================================================
# GAMMA DIVERSITY CHANGE (1000m, 2024-2023)
# ============================================================================

gamma_change <- master_gamma_results_2024 %>%
  dplyr::group_by(Transect, SoilVeg, Treatment) %>%
  dplyr::mutate(replicate = row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::rename(sr_2024 = sr) %>%
  left_join(
    master_gamma_results_2023 %>%
      dplyr::group_by(Transect, SoilVeg, Treatment) %>%
      dplyr::mutate(replicate = row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::select(Transect, SoilVeg, Treatment, replicate, sr_2023 = sr),
    by = c("Transect", "SoilVeg", "Treatment", "replicate")
  ) %>%
  dplyr::mutate(sr_change = sr_2024 - sr_2023) %>%
  na.omit()

ggplot(gamma_change, aes(Treatment, sr_change, fill = Treatment))+
  facet_wrap(~SoilVeg)+
  geom_boxplot()+
  scale_fill_manual(values = c( "#FFE793", "#4F93A7"))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  ylim(-20, 20)+
  xlab("")+
  ylab("Change in gamma diversity (1,000 m2 richness)")+
  theme_base()

ggsave("C:/Users/ohler/Dropbox/grants/Gemini/figures/gammadiv_change.pdf",
       plot = last_plot(),
       device = "pdf",
       path = NULL,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in"),
       dpi = 600,
       limitsize = TRUE)

mod <- lme(sr_change~Treatment*SoilVeg, random = ~1|Transect, data = gamma_change)
summary(mod)
emmeans(mod, ~ Treatment*SoilVeg)
pairs(emmeans(mod, ~ Treatment*SoilVeg))

