####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 04-removal-analysis.R
# Created April 2022
# Last Updated November 2022

####### Import Libraries and External Files #######

library(cmdstanr)
library(napops)
library(plyr)
library(ggplot2)
library(ggpubr)
theme_set(theme_pubclean())

####### Read Data #################################

brownian <- readRDS("output/cv/removal/brownian/1full.RDS")
load("data/generated/removal_stan_data_cv.rda")
load("data/generated/corr_matrix_cv.rda")
binomial <- read.csv("data/generated/binomial_names.csv")

####### Main Code #################################

# Extract log_phi summary statistics from full Stan model runs
brownian_summary <- brownian$summary(variables = "log_phi")
# pagel_summary <- pagel$summary(variables = "log_phi")

# Add species names to these summaries
brownian_summary$Scientific_BT <- gsub("_", " ", rownames(corr_matrix_cv))
brownian_summary <- join(brownian_summary, binomial[, c("Scientific_BT", "Code")], by = "Scientific_BT")

# pagel_summary$Scientific_BT <- gsub("_", " ", rownames(corr_matrix_cv))
# pagel_summary <- join(pagel_summary, binomial[, c("Scientific_BT", "Code")], by = "Scientific_BT")

# Get original single-species NA-POPS estimates
napops_summary <- napops::coef_removal(species = brownian_summary$Code, model = 1)

to_plot <- merge(brownian_summary[,c("Code", "mean")], napops_summary[, c("Species", "Intercept")],
                 by.x = "Code", by.y = "Species")
names(to_plot) <- c("Species", "Multi", "Single")

ggplot(data = to_plot, mapping = aes(x = exp(Single), y = exp(Multi))) +
  geom_point() +
  geom_abline(slope = 1) +
  xlim(0,1) +
  ylim(0,1) +
  NULL

to_plot$CR_Multi <- exp(to_plot$Multi)
to_plot$CR_Single <- exp(to_plot$Single)
to_plot$difference <- to_plot$CR_Multi - to_plot$CR_Single

to_plot <- to_plot[match(brownian_summary$Code, to_plot$Species), ]
to_plot$Species <- factor(to_plot$Species, levels = c(to_plot$Species))

png(filename = "output/plots/removal_diff_brownian.png",
    width = 6, height = 40, res = 300, units = "in")
ggplot(data = to_plot, mapping = aes(x = Species, y = difference)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  NULL
dev.off()

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
