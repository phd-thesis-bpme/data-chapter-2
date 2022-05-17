####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 04-removal-analysis.R
# Created April 2022
# Last Updated April 2022

####### Import Libraries and External Files #######

library(rstan)
library(napops)
library(ggplot2)
library(ggpubr)
theme_set(theme_pubclean())

####### Read Data #################################

load("output/turdidae_corr_matrix.rda")
load("output/removal_turdidae_model.rda")
load("output/removal_turdidae_data.rda")
binomial <- read.csv("output/binomial_names.csv")

####### Main Code #################################

species <- rownames(corr_matrix)
species <- gsub("_", " ", species)
bin_red <- binomial[which(binomial$Scientific %in% species), ]
bin_red <- bin_red[match(species, bin_red$Scientific), ]
codes <- bin_red$Code

# Create empty data frame
removal_df <- data.frame(Species = rep(codes, 2),
                         N = NA,
                         Estimate = NA,
                         Lower = NA,
                         Upper = NA,
                         Model = c(rep("Single", 11), rep("Multi", 11)))
stan_summary <- data.frame(summary(stan_job)$summary)

i <- 1
for (sp in codes)
{
  napops_est <- napops::cue_rate(species = sp, model = 1, od = c(100:105), tssr = c(-1:1), quantiles = c(0.025, 0.975))
  removal_df[which(removal_df$Species == sp & 
                     removal_df$Model == "Single"), "Estimate"] <- napops_est$CR_est[1]
  removal_df[which(removal_df$Species == sp & 
                     removal_df$Model == "Single"), "Lower"] <- napops_est$CR_2.5[1]
  removal_df[which(removal_df$Species == sp & 
                     removal_df$Model == "Single"), "Upper"] <- napops_est$CR_97.5[1]
  removal_df[which(removal_df$Species == sp & 
                     removal_df$Model == "Single"), "N"] <- nrow(covariates_removal(species = sp))
  removal_df[which(removal_df$Species == sp & 
                     removal_df$Model == "Multi"), "N"] <- nrow(covariates_removal(species = sp))
  
  removal_df[which(removal_df$Species == sp & 
                     removal_df$Model == "Multi"), "Estimate"] <- exp(stan_summary$mean[i])
  removal_df[which(removal_df$Species == sp & 
                     removal_df$Model == "Multi"), "Lower"] <- exp(stan_summary$X2.5.[i])
  removal_df[which(removal_df$Species == sp & 
                     removal_df$Model == "Multi"), "Upper"] <- exp(stan_summary$X97.5.[i])
  
  i <- i + 1
}

removal_df$Species <- factor(removal_df$Species, levels = codes)
removal_df$Model <- factor(removal_df$Model, levels = c("Single", "Multi"))
removal_df$Log_N <- log(removal_df$N)

# Plot single vs multi
comparison_plot <- ggplot(data = removal_df, aes(x = Species, y = Estimate, group = Model, color = Model)) +
  geom_point(aes(size = Log_N), position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, position = position_dodge(width = 0.4)) +
  ylim(0,0.5) +
  ylab("Cue Rate") +
  NULL


####### Output ####################################

png("plots/removal_comparison.png",
    width = 6, height = 4, res = 600, units = "in")
print(comparison_plot)
dev.off()
