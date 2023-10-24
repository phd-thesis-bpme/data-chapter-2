####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 8-removal-comparison-analysis.R
# Created April 2022
# Last Updated October 2023

####### Import Libraries and External Files #######

library(cmdstanr)
library(bayesplot)
library(napops)
library(plyr)
library(ggplot2)
library(ggpubr)
theme_set(theme_pubclean())
bayesplot::color_scheme_set("red")

####### Read Data #################################

rem_model <- readRDS("output/model_runs/removal_predictions.RDS")
load("data/generated/removal_stan_data_pred.rda")

####### Data Wrangling ############################

# Extract log_phi summary statistics from full Stan model runs
rem_summary <- rem_model$summary("log_phi")

# Add species names to these summaries
rem_summary$Code <- removal_stan_data$sp_all

# Get data sample size for all species and add to summary
species_n <- data.frame(table(removal_stan_data$species))
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
  geom_text(aes(label = Label)) +
  xlab("Cue Rate (Single Species)") + 
    ylab("Cue Rate (Multi Species)") +
  NULL)

(difference_plot <- ggplot(data = to_plot, mapping = aes(x = Species, y = sort(diff))) +
    geom_point() + 
    geom_text(aes(label = Label)) +
    geom_hline(yintercept = 0) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ylab("Difference in Cue Rate (Multi - Single Model)") +
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

####### SD Comparison Plot ########################

to_plot <- rem_summary[which(rem_summary$Code %in% napops_summary$Species), 
                       c("Code", "N", "sd")]
names(to_plot) <- c("Code", "N", "SD_Multi")
to_plot$SD_Single <- NA
# Get original NA-POPS SDs
for (i in 1:nrow(to_plot))
{
  to_plot$SD_Single[i] <- sqrt(napops::get_vcv(species = to_plot$Code[i],
                                          model_type = "rem",
                                          model_num = 1))
}
to_plot_long <- reshape2::melt(to_plot, id = c("Code", "N"), variable_name = "Model")

sd_model_data <- list(n_samples = nrow(to_plot_long),
                      stdev = to_plot_long$value,
                      N = to_plot_long$N,
                      model_type = factor(to_plot_long$variable,
                                     levels = c("SD_Single", "SD_Multi")))

sd_model <- cmdstan_model(stan_file = "models/sd_model.stan")

sd_model_run <- sd_model$sample(
  data = sd_model_data,
  iter_warmup = 1000,
  iter_sampling = 2000,
  chains = 4,
  parallel_chains = 4,
  refresh = 10
)

(sd_intercept_plot <- bayesplot::mcmc_areas(sd_model_run$draws(c("alpha[1]", "alpha[2]")),
                      prob = 0.95))

(sd_slope_plot <- bayesplot::mcmc_areas(sd_model_run$draws(c("beta[1]", "beta[2]")),
                      prob = 0.95))

####### Species-specific Plots ####################

sp <- c("LCTH", "LEPC", "HASP", "TRBL", "SPOW", "KIWA", "BITH")

to_plot <- rem_summary[which(rem_summary$Code %in% sp), ]

species_vars <- to_plot$variable

(species_prediction_plot <- bayesplot::mcmc_intervals(exp(rem_model$draws(species_vars))) + 
  xlab("Predicted Cue Rate") + 
  ylab("Species") +
  scale_y_discrete(labels = (to_plot$Code)) +
  coord_flip())

####### Plotting Other Parameters #################

(intercept_plot <- bayesplot::mcmc_areas(rem_model$draws(c("intercept")),
                                         prob = 0.95))

(mig_strat_plot <- bayesplot::mcmc_areas(rem_model$draws(c("mu_mig_strat[1]", "mu_mig_strat[2]")),
                                        prob = 0.95))

(sigma_plot <- bayesplot::mcmc_areas(rem_model$draws(c("sigma")),
                                         prob = 0.95))

####### Output ####################################

png("output/plots/removal_1vs1.png",
    width = 6, height = 3, res = 600, units = "in")
ggarrange(modelled_difference_plot, single_vs_multi_plot, ncol = 2, labels = c("A", "B"))
dev.off()

png("output/plots/removal_sd_plot.png",
    width = 6, height = 4, res = 600, units = "in")
ggarrange(sd_intercept_plot, sd_slope_plot, ncol = 2, labels = c("A", "B"))
dev.off()

png("output/plots/removal_predictions_plot.png",
    width = 6, height = 4, res = 600, units = "in")
print(species_prediction_plot)
dev.off()

png("output/plots/removal_other_params.png",
    width = 6, height = 6, res = 600, units = "in")
ggarrange(intercept_plot, mig_strat_plot, sigma_plot, nrow = 3, labels = c("A", "B", "C"))
dev.off()

