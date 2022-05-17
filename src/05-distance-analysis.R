####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 05-distance-analysis.R
# Created May 2022
# Last Updated May 2022

####### Import Libraries and External Files #######

library(rstan)
library(napops)
library(ggplot2)
library(ggpubr)
theme_set(theme_pubclean())

####### Read Data #################################

load("output/turdidae_corr_matrix.rda")
load("output/distance_turdidae_model.rda")
load("output/distance_turdidae_data.rda")
binomial <- read.csv("output/binomial_names.csv")

####### Main Code #################################

species <- rownames(corr_matrix)
species <- gsub("_", " ", species)
bin_red <- binomial[which(binomial$Scientific %in% species), ]
bin_red <- bin_red[match(species, bin_red$Scientific), ]
codes <- bin_red$Code

# Create empty data frame
distance_df <- data.frame(Species = rep(codes, 2),
                         N = NA,
                         Estimate = NA,
                         Lower = NA,
                         Upper = NA,
                         Model = c(rep("Single", 11), rep("Multi", 11)))
stan_summary <- data.frame(summary(stan_job)$summary)

i <- 1
for (sp in codes)
{
  napops_est <- napops::edr(species = sp, model = 1, forest = c(0.0,0.5,1.0), road = TRUE, quantiles = c(0.025, 0.975))
  distance_df[which(distance_df$Species == sp & 
                     distance_df$Model == "Single"), "Estimate"] <- napops_est$EDR_est[1]
  distance_df[which(distance_df$Species == sp & 
                     distance_df$Model == "Single"), "Lower"] <- napops_est$EDR_2.5[1]
  distance_df[which(distance_df$Species == sp & 
                     distance_df$Model == "Single"), "Upper"] <- napops_est$EDR_97.5[1]
  distance_df[which(distance_df$Species == sp & 
                     distance_df$Model == "Single"), "N"] <- nrow(covariates_distance(species = sp))
  distance_df[which(distance_df$Species == sp & 
                     distance_df$Model == "Multi"), "N"] <- nrow(covariates_distance(species = sp))
  
  distance_df[which(distance_df$Species == sp & 
                     distance_df$Model == "Multi"), "Estimate"] <- exp(stan_summary$mean[i])
  distance_df[which(distance_df$Species == sp & 
                     distance_df$Model == "Multi"), "Lower"] <- exp(stan_summary$X2.5.[i])
  distance_df[which(distance_df$Species == sp & 
                     distance_df$Model == "Multi"), "Upper"] <- exp(stan_summary$X97.5.[i])
  
  i <- i + 1
}

distance_df$Species <- factor(distance_df$Species, levels = codes)
distance_df$Model <- factor(distance_df$Model, levels = c("Single", "Multi"))
distance_df$Log_N <- log(distance_df$N)

# Plot single vs multi
comparison_plot <- ggplot(data = distance_df, aes(x = Species, y = Estimate, group = Model, color = Model)) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, position = position_dodge(width = 0.4)) +
  ylab("EDR") +
  NULL


####### Output ####################################

png("plots/distance_comparison.png",
    width = 6, height = 4, res = 600, units = "in")
print(comparison_plot)
dev.off()
