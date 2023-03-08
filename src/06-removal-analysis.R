####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 06-removal-analysis.R
# Created April 2022
# Last Updated January 2023

####### Import Libraries and External Files #######

library(cmdstanr)
library(bayesplot)
library(napops)
library(plyr)
library(ggplot2)

####### Read Data #################################

rem_model <- readRDS("output/model_runs/removal_predictions.RDS")
load("data/generated/corr_matrix_predict.rda")
binomial <- read.csv("data/generated/binomial_names.csv")
load("data/generated/removal_stan_data_pred.rda")

####### Data Wrangling ############################

# Extract log_phi summary statistics from full Stan model runs
rem_summary <- rem_model$summary(variables = "log_phi")

# Model diagnostics

diag <- rem_model$sampler_diagnostics(format = "df")
sigma_draws <- rem_model$draws(variables = "mu_mig_strat[1]", format = "df")
sigma_draws$.chain <- factor(sigma_draws$.chain, levels = c("1", "2", "3", "4"))
ggplot(data = sigma_draws, aes(x = .iteration, y = mu, group = (.chain), color = .chain)) +
  geom_line()
mcmc_hist(rem_model$draws("sigma")) 

# Add species names to these summaries
rem_summary$Scientific_BT <- gsub("_", " ", rownames(corr_matrix_predict))
rem_summary <- join(rem_summary, binomial[, c("Scientific_BT", "Code")], by = "Scientific_BT")

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

to_plot <- merge(rem_summary[,c("Code", "mean")], napops_summary[, c("Species", "Intercept")],
                 by.x = "Code", by.y = "Species")
names(to_plot) <- c("Species", "Multi", "Single")
to_plot$Label <- ""
to_plot$diff <- NA
for (i in 1:nrow(to_plot))
{
  to_plot$diff[i] <- abs(to_plot$Multi[i] - to_plot$Single[i])
  if (to_plot$diff[i] > 0.1)
  {
    to_plot$Label[i] <- to_plot$Species[i]
  }
}
single_vs_multi_plot <- ggplot(data = to_plot, mapping = aes(x = exp(Single), y = exp(Multi))) +
  geom_point() +
  geom_abline(slope = 1) +
  xlim(0,1) +
  ylim(0,1) +
  geom_text(aes(label = Label)) +
  NULL

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

sd_comp_plot <- ggplot(data = to_plot_long, 
                       aes(x = log(N), y = log(value),
                           group = variable, color = variable)) + 
  geom_smooth() +
  NULL

####### Species-specific Plots ####################

sp <- c("LCTH", "LEPC", "HASP", "TRBL", "SPOW", "KIWA", "BITH")

to_plot <- rem_summary[which(rem_summary$Code %in% sp), ]

species_plot <- ggplot(data = to_plot, aes(x = Code, y = exp(mean),
                                           ymin = exp(q5), ymax = exp(q95))) +
  geom_point() + 
  geom_errorbar() +
  xlab("Species") +
  ylab("Cue Rate") +
  geom_text(aes(y = exp(q95) + 0.25, label = N)) +
  NULL

####### Output ####################################

png("output/plots/removal_1vs1.png",
    width = 6, height = 4, res = 600, units = "in")
print(single_vs_multi_plot)
dev.off()

png("output/plots/removal_sd_plot.png",
    width = 6, height = 4, res = 600, units = "in")
print(sd_comp_plot)
dev.off()

png("output/plots/removal_predictions_plot.png",
    width = 6, height = 4, res = 600, units = "in")
print(species_plot)
dev.off()
