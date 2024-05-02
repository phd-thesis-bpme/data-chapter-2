####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# posthoc/2-sd-figure.R
# Created December 2023
# Last Updated April 2024

####### Import Libraries and External Files #######

library(cmdstanr)
library(bayesplot)
library(plyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(patchwork)
theme_set(theme_pubclean())
bayesplot::color_scheme_set("red")

####### Read Data #################################

rem_ms <- readRDS("output/model_runs/removal_ms.RDS")
rem_ss <- readRDS("output/model_runs/removal_ss.RDS")
load("data/generated/removal_stan_data_cv.rda")

dis_ms <- readRDS("output/model_runs/distance_ms.RDS")
dis_ss <- readRDS("output/model_runs/distance_ss.RDS")
load("data/generated/distance_stan_data_cv.rda")

####### Removal Model #############################

# Extract log_phi summary statistics from full Stan model runs
rem_summary_ms <- rem_ms$summary("log_phi")
rem_summary_ss <- rem_ss$summary("log_phi")

# Add species names to these summaries
rem_summary_ms$Code <- removal_stan_data_cv$sp_all
rem_summary_ss$Code <- removal_stan_data_cv$sp_all

# Get data sample size for all species and add to summary
species_n <- data.frame(table(removal_stan_data_cv$species))
names(species_n) <- c("index", "N")
rem_summary_ms$index <- seq(1, nrow(rem_summary_ms))
rem_summary_ss$index <- seq(1, nrow(rem_summary_ss))

rem_summary_ms$N <- 0
rem_summary_ss$N <- 0

for (i in 1:nrow(species_n))
{
  rem_summary_ms[which(rem_summary_ms$index == species_n$index[i]), "N"] <-
    species_n$N[i]
  
  rem_summary_ss[which(rem_summary_ss$index == species_n$index[i]), "N"] <-
    species_n$N[i]
}

to_plot <- data.frame(Species = rep(removal_stan_data_cv$sp_all, 2),
                      N = rep(species_n$N, 2),
                      SD = c(rem_summary_ms$sd, rem_summary_ss$sd),
                      Model = c(rep("Multi", nrow(rem_summary_ms)),
                                rep("Single", nrow(rem_summary_ss))))

(sd_comp_removal <- ggplot(data = to_plot, 
                       aes(x = log(N), y = SD,
                           group = Model, color = Model)) + 
  geom_smooth() +
  xlab("log(Sample Size)") +
  ylab("Standard Deviation") +
  theme(legend.position = "none") +
  NULL)

####### Distance Model #############################

# Extract log_tau summary statistics from full Stan model runs
dis_summary_ms <- dis_ms$summary("log_tau")
dis_summary_ss <- dis_ss$summary("log_tau")

# Add the species names
dis_summary_ms$Code <- distance_stan_data_cv$sp_all
dis_summary_ss$Code <- distance_stan_data_cv$sp_all

# Get data sample size for all species and add to summary
species_n <- data.frame(table(distance_stan_data_cv$species))
names(species_n) <- c("index", "N")
dis_summary_ms$index <- seq(1, nrow(dis_summary_ms))
dis_summary_ss$index <- seq(1, nrow(dis_summary_ss))

dis_summary_ms$N <- 0
dis_summary_ss$N <- 0

for (i in 1:nrow(species_n))
{
  dis_summary_ms[which(dis_summary_ms$index == species_n$index[i]), "N"] <-
    species_n$N[i]
  
  dis_summary_ss[which(dis_summary_ss$index == species_n$index[i]), "N"] <-
    species_n$N[i]
}

to_plot <- data.frame(Species = rep(distance_stan_data_cv$sp_all, 2),
                      N = rep(species_n$N, 2),
                      SD = c(dis_summary_ms$sd, dis_summary_ss$sd),
                      Model = c(rep("Multi", nrow(dis_summary_ms)),
                                rep("Single", nrow(dis_summary_ss))))

(sd_comp_distance <- ggplot(data = to_plot, 
                           aes(x = log(N), y = SD,
                               group = Model, color = Model)) + 
    geom_smooth() +
    xlab("log(Sample Size)") +
    ylab("Standard Deviation") +
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

