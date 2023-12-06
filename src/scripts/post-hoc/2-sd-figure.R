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

n_vals <- seq(min(to_plot_long$N), max(to_plot_long$N))

ss_draws <- as.matrix(sd_model_run$draws(c("alpha[1]", "beta[1]"), format = "df"))[, 1:2]
ss_pred <- matrix(data = NA, nrow = length(n_vals), ncol = nrow(ss_draws))
for (i in 1:nrow(ss_draws))
{
  ss_pred[ ,i] <- ss_draws[i, 1] + ss_draws[i, 2] * log(n_vals)
}
ss_pred_mean <- apply(ss_pred, 1, FUN = mean)
ss_q25 <- apply(ss_pred, 1, FUN = quantile, probs = 0.025)
ss_q975 <- apply(ss_pred, 1, FUN = quantile, probs = 0.975)
rm(ss_pred); gc()

ms_draws <- as.matrix(sd_model_run$draws(c("alpha[2]", "beta[2]"), format = "df"))[, 1:2]
ms_pred <- matrix(data = NA, nrow = length(n_vals), ncol = nrow(ms_draws))
for (i in 1:nrow(ms_draws))
{
  ms_pred[ ,i] <- ms_draws[i, 1] + ms_draws[i, 2] * log(n_vals)
}
ms_pred_mean <- apply(ms_pred, 1, FUN = mean)
ms_q25 <- apply(ms_pred, 1, FUN = quantile, probs = 0.025)
ms_q975 <- apply(ms_pred, 1, FUN = quantile, probs = 0.975)
rm(ms_pred); gc()

to_plot_preds <- data.frame(N = rep(n_vals, 2),
                      Mean = c(ss_pred_mean, ms_pred_mean),
                      q25 = c(ss_q25, ms_q25),
                      q975 = c(ss_q975, ms_q975),
                      Model = c(rep("Single", length(n_vals)),
                                rep("Multi", length(n_vals))))

sd_comp_plot <- ggplot(to_plot_preds, aes(x = log(N), y = Mean,
                                          ymin = q25, ymax = q975,
                                          color = as.factor(Model))) +
  geom_line() + 
  geom_ribbon(alpha = 0.25)

(sd_intercept_plot <- bayesplot::mcmc_areas(sd_model_run$draws(c("alpha[1]", "alpha[2]")),
                                            prob = 0.95))

(sd_slope_plot <- bayesplot::mcmc_areas(sd_model_run$draws(c("beta[1]", "beta[2]")),
                                        prob = 0.95))