####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# posthoc/2-sd-figure.R
# Created December 2023
# Last Updated December 2023

####### Import Libraries and External Files #######

library(cmdstanr)
library(bayesplot)
library(napops)
library(plyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(patchwork)
theme_set(theme_pubclean())
bayesplot::color_scheme_set("red")

####### Read Data #################################

rem_model <- readRDS("output/model_runs/removal_predictions.RDS")
load("data/generated/removal_stan_data_pred.rda")
dis_model <- readRDS("output/model_runs/distance_predictions.RDS")
load("data/generated/distance_stan_data.rda")

####### Removal Model #############################

# Extract log_phi summary statistics from full Stan model runs
rem_summary <- rem_model$summary("log_phi")

# Add species names to these summaries
rem_summary$Code <- removal_stan_data_pred$sp_all

# Get data sample size for all species and add to summary
species_n <- data.frame(table(removal_stan_data_pred$species))
names(species_n) <- c("index", "N")
rem_summary$index <- seq(1, nrow(rem_summary))
rem_summary$N <- 0
for (i in 1:nrow(species_n))
{
  rem_summary[which(rem_summary$index == species_n$index[i]), "N"] <-
    species_n$N[i]
}

# Get original single-species NA-POPS estimates
napops_summary <- napops::coef_removal(species = rem_summary$Code, model = 1)

to_plot <- rem_summary[which(rem_summary$Code %in% napops_summary$Species), 
                       c("Code", "N", "sd")]
names(to_plot) <- c("Code", "N", "Multi")
to_plot$SD_Single <- NA
# Get original NA-POPS SDs
for (i in 1:nrow(to_plot))
{
  to_plot$SD_Single[i] <- sqrt(napops::get_vcv(species = to_plot$Code[i],
                                               model_type = "rem",
                                               model_num = 1))
}
to_plot_long <- reshape2::melt(to_plot, id = c("Code", "N"), variable_name = "Model")
(sd_comp_removal <- ggplot(data = to_plot_long, 
                       aes(x = log(N), y = (value),
                           group = variable, color = variable)) + 
  geom_smooth() +
  xlab("log(Sample Size)") +
  ylab("Standard Deviation") +
    theme(legend.position = "none") +
  NULL)

####### Distance Model #############################

# Extract log_phi summary statistics from full Stan model runs
dis_summary <- dis_model$summary("log_tau")

# Add species names to these summaries
dis_summary$Code <- distance_stan_data$sp_all

# Get data sample size for all species and add to summary
species_n <- data.frame(table(distance_stan_data$species))
names(species_n) <- c("index", "N")
dis_summary$index <- seq(1, nrow(dis_summary))
dis_summary$N <- 0
for (i in 1:nrow(species_n))
{
  dis_summary[which(dis_summary$index == species_n$index[i]), "N"] <-
    species_n$N[i]
}

# Get original single-species NA-POPS estimates
napops_summary <- napops::coef_distance(species = dis_summary$Code, model = 1)

to_plot <- dis_summary[which(dis_summary$Code %in% napops_summary$Species), 
                       c("Code", "N", "sd")]
names(to_plot) <- c("Code", "N", "SD_Multi")
to_plot$SD_Single <- NA
# Get original NA-POPS SDs
for (i in 1:nrow(to_plot))
{
  to_plot$SD_Single[i] <- sqrt(napops::get_vcv(species = to_plot$Code[i],
                                               model_type = "dis",
                                               model_num = 1))
}
to_plot_long <- reshape2::melt(to_plot, id = c("Code", "N"), variable_name = "Model")
(sd_comp_distance <- ggplot(data = to_plot_long, 
                           aes(x = log(N), y = (value),
                               group = variable, color = variable)) + 
    geom_smooth() +
    xlab("log(Sample Size)") +
    ylab("Standard Deviation") +
    ylim(0,0.4) +
    theme(legend.position = "none") +
    NULL)

####### Output ####################################

tiff("output/plots/sd_comp.tiff",
     width = 6, height = 3, res = 600, units = "in")
ggarrange(sd_comp_removal, sd_comp_distance, ncol = 2, labels = c("A", "B"))
dev.off()

png("output/plots/sd_comp.png",
     width = 6, height = 3, res = 600, units = "in")
ggarrange(sd_comp_removal, sd_comp_distance, ncol = 2, labels = c("A", "B"))
dev.off()

