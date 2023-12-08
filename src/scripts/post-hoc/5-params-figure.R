####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# posthoc/5-params-figure.R
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
dis_model <- readRDS("output/model_runs/distance_predictions.RDS")

####### Removal Model #############################

(removal_plot <- bayesplot::mcmc_intervals(rem_model$draws(c("sigma",
                                                             "mu_mig_strat[2]", 
                                                             "mu_mig_strat[1]",
                                                             "intercept"))) +
    scale_y_discrete(labels = c("Sigma", "Migrant", "Resident", "Intercept")))

params_removal <- rem_model$summary(c("intercept", "mu_mig_strat", "sigma"))

####### Distance Model ############################

(distance_plot <- bayesplot::mcmc_intervals(dis_model$draws(c("sigma",
                                                              "beta_pitch",
                                                              "beta_mass", 
                                                              "mu_habitat[2]",
                                                              "mu_habitat[1]",
                                                              "mu_mig_strat[2]",
                                                              "mu_mig_strat[1]",
                                                              "intercept"))) +
   scale_y_discrete(labels = c("Sigma", "Pitch Slope","Mass Slope","Closed Habitat",
                               "Open Habitat","Migrant","Resident","Intercept")))

params_distance <- dis_model$summary(c("intercept", "mu_mig_strat",
                                      "mu_habitat", "beta_mass", "beta_pitch",
                                      "sigma"))

####### Output ####################################

write.table(params_removal, file = "data/generated/removal_params.csv", sep = ",", row.names = FALSE)
write.table(params_distance, file = "data/generated/distance_params.csv", sep = ",", row.names = FALSE)

tiff(filename = "output/plots/parameters_plot.tiff",
     width = 6, height = 6, units = "in", res = 600)
ggarrange(removal_plot, distance_plot, nrow = 2, labels = c("A", "B"))
dev.off()

png(filename = "output/plots/parameters_plot.png",
     width = 6, height = 6, units = "in", res = 600)
ggarrange(removal_plot, distance_plot, nrow = 2, labels = c("A", "B"))
dev.off()
