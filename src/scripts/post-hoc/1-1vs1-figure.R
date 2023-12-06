####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# posthoc/1-1vs1-figure.R
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

####### 1-to-1 Plot ###############################

to_plot <- merge(rem_summary[,c("Code", "mean", "N")], napops_summary[, c("Species", "Intercept")],
                 by.x = "Code", by.y = "Species")
names(to_plot) <- c("Species", "Multi", "N", "Single")
to_plot$Label <- ""
to_plot$Multi <- exp(to_plot$Multi)
to_plot$Single <- exp(to_plot$Single)
to_plot$diff <- to_plot$Multi - to_plot$Single

# Check for a 20% relative difference in cue rates
for (i in 1:nrow(to_plot))
{
  if ((abs(to_plot$diff[i]) / to_plot$Single[i]) > 0.20)
  {
    to_plot$Label[i] <- to_plot$Species[i]
  }
}

(single_vs_multi_plot <- ggplot(data = to_plot, mapping = aes(x = Single, y = Multi)) +
    geom_point() +
    geom_abline(slope = 1) +
    xlim(0,1) +
    ylim(0,1) +
    geom_text_repel(box.padding = 0.5, max.overlaps = Inf, aes(label = Label)) +
    xlab("Cue Rate (Single Species)") + 
    ylab("Cue Rate (Multi Species)") +
    NULL)

diff_model_data <- list(n_samples = nrow(to_plot),
                        difference = to_plot$diff,
                        N = to_plot$N)
diff_model <- cmdstan_model(stan_file = "models/difference.stan")

diff_model_run <- diff_model$sample(
  data = diff_model_data,
  iter_warmup = 1000,
  iter_sampling = 2000,
  chains = 4,
  parallel_chains = 4,
  refresh = 10
)

(modelled_difference_plot <- bayesplot::mcmc_areas(diff_model_run$draws(c("mu")),
                                                   prob = 0.95) +
    xlab("Modelled Difference"))

removal_plot <- single_vs_multi_plot + inset_element(modelled_difference_plot,
                                                     left = 0.6, 
                                                     bottom = 0, 
                                                     right = 1, 
                                                     top = 0.5)

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

to_plot <- merge(dis_summary[,c("Code", "mean")], napops_summary[, c("Species", "Intercept")],
                 by.x = "Code", by.y = "Species")
names(to_plot) <- c("Species", "Multi", "Single")
to_plot$Label <- ""
to_plot$Multi <- exp(to_plot$Multi) * 100
to_plot$Single <- exp(to_plot$Single)
to_plot$diff <- to_plot$Multi - to_plot$Single

# Check for a 20% relative difference in EDR
for (i in 1:nrow(to_plot))
{
  if ((abs(to_plot$diff[i]) / to_plot$Single[i]) > 0.20)
  {
    to_plot$Label[i] <- to_plot$Species[i]
  }
}

(single_vs_multi_plot <- ggplot(data = to_plot, mapping = aes(x = Single, y = Multi)) +
    geom_point() +
    geom_abline(slope = 1) +
    xlim(0,500) +
    ylim(0,500) +
    geom_text_repel(box.padding = 0.5, max.overlaps = Inf, aes(label = Label)) +
    xlab("EDR (Single Species)") + 
    ylab("EDR (Multi Species)") +
    NULL)

diff_model_data <- list(n_samples = nrow(to_plot),
                        difference = to_plot$diff)
diff_model <- cmdstan_model(stan_file = "models/difference.stan")

diff_model_run <- diff_model$sample(
  data = diff_model_data,
  iter_warmup = 1000,
  iter_sampling = 2000,
  chains = 4,
  parallel_chains = 4,
  refresh = 10
)

(modelled_difference_plot <- bayesplot::mcmc_areas(diff_model_run$draws(c("mu")),
                                                   prob = 0.95) +
    xlab("Modelled Difference"))

distance_plot <- single_vs_multi_plot + inset_element(modelled_difference_plot,
                                                     left = 0.6, 
                                                     bottom = 0, 
                                                     right = 1, 
                                                     top = 0.5)

####### Output ####################################

tiff("output/plots/1vs1.tiff",
     width = 6, height = 6, res = 600, units = "in")
ggarrange(removal_plot, distance_plot, nrow = 2, labels = c("A", "B"))
dev.off()
