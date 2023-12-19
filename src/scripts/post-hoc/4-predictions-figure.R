####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# posthoc/4-predictions-figure.R
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
library(gginnards)
library(patchwork)
theme_set(theme_pubclean())
bayesplot::color_scheme_set("red")

####### Read Data #################################

rem_model <- readRDS("output/model_runs/removal_predictions.RDS")
load("data/generated/removal_stan_data_pred.rda")
dis_model <- readRDS("output/model_runs/distance_predictions.RDS")
load("data/generated/distance_stan_data.rda")

traits <- read.csv("data/raw/traits.csv")
binomial <- read.csv("data/generated/binomial_names.csv")

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
rem_summary <- dplyr::left_join(x = rem_summary, y = binomial[, c("Code", "Scientific_BT")],
                                by = "Code")
rem_summary$Scientific_BT <- gsub(x = rem_summary$Scientific_BT, pattern = " ", replacement = "_")

rem_summary <- dplyr::left_join(x = rem_summary, y = traits[, c("Code", "Migrant")],
                                 by = "Code")

sp <- c("LCTH", "LEPC", "HASP", "TRBL", "SPOW", "KIWA", "BITH")

to_plot <- rem_summary[which(rem_summary$Code %in% sp), ]

species_vars <- to_plot$variable

other_sp_df <- data.frame(Species = character(),
                          mean = numeric(),
                          Alpha = numeric())
phylo <- removal_stan_data_pred$phylo_corr
for (sv in species_vars)
{
  sci_name <- rem_summary[which(rem_summary$variable == sv), ]$Scientific_BT
  migrant <- rem_summary[which(rem_summary$variable == sv), ]$Migrant
  
  alpha_df <- data.frame(Sci_Name = dimnames(phylo)[[1]],
                         Alpha = phylo[sci_name, ])
  

  alpha_df <- merge(alpha_df, rem_summary[which(rem_summary$Migrant == migrant &
                                                  (rem_summary$Code %in% sp) == FALSE), 
                                          c("mean", "Scientific_BT")],
                    by.x = "Sci_Name", by.y = "Scientific_BT")
  alpha_df$mean <- exp(alpha_df$mean)
  alpha_df$Species <- rep(sv, times = nrow(alpha_df))
  
  other_sp_df <- rbind(other_sp_df,
                       alpha_df[, c("Species", "mean", "Alpha")])
}

other_sp_df <- other_sp_df[which(other_sp_df$Alpha > 0.001), ]

rem_draws <- exp(rem_model$draws(species_vars))
attributes(rem_draws)$dimnames$variable <- to_plot$Code

(removal_plot <- bayesplot::mcmc_intervals(exp(rem_model$draws(species_vars))) + 
    xlab("Predicted Cue Rate") + 
    ylab("Species") +
    #geom_point(data = test_df, aes(x = exp(cr), y = Species,  alpha = alpha), position= "jitter") +
    #scale_y_discrete(labels = (to_plot$Code)) +
    coord_flip())

sp_ordered <- rev(to_plot$Code)
removal_plot$layers <- c(geom_point(data = other_sp_df, aes(x = mean, y = Species,  alpha = Alpha), position=position_jitter(height=0.2)),
                         removal_plot$layers)
removal_plot <- removal_plot + scale_y_discrete(labels = rev(to_plot$Code))


removal_plot
                

removal_plot <- removal_plot + scale_y_discrete(labels = (to_plot$Code))


# Extract mean cue rates of "similar" species
sp_similar <- binomial[which(binomial$Family %in% c("Passerellidae",
                                                    "Strigidae",
                                                    "Turdidae",
                                                    "Mimidae",
                                                    "Icteridae",
                                                    "Parulidae",
                                                    "Phasianidae"))]



####### Distance Model ############################

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

sp <- c("LCTH", "LEPC", "HASP", "SPOW", "KIWA", "BITH")

to_plot <- dis_summary[which(dis_summary$Code %in% sp), ]

species_vars <- to_plot$variable

(distance_plot <- bayesplot::mcmc_intervals(exp(dis_model$draws(species_vars)) * 100) + 
    xlab("Predicted EDR") + 
    ylab("Species") +
    scale_y_discrete(labels = (to_plot$Code)) +
    coord_flip())

####### Output ####################################

tiff(filename = "output/plots/predictions_figure.tiff",
     width = 6, height = 6, units = "in", res = 600)
ggarrange(removal_plot, distance_plot, nrow = 2, labels = c("A", "B"))
dev.off()

png(filename = "output/plots/predictions_figure.png",
     width = 6, height = 6, units = "in", res = 600)
ggarrange(removal_plot, distance_plot, nrow = 2, labels = c("A", "B"))
dev.off()
