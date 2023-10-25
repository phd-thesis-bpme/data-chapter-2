####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 8-distance-comparison-analysis.R
# Created May 2022
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

dis_model <- readRDS("output/model_runs/distance_predictions.RDS")
load("data/generated/distance_stan_data.rda")

####### Data Wrangling ############################

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

####### 1-to-1 Plot ###############################

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
    geom_text(aes(label = Label)) +
    xlab("EDR (Single Species)") + 
    ylab("EDR (Multi Species)") +
    NULL)

(difference_plot <- ggplot(data = to_plot, mapping = aes(x = Species, y = sort(diff))) +
    geom_point() + 
    geom_text(aes(label = Label)) +
    geom_hline(yintercept = 0) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ylab("Difference in EDR (Multi - Single Model)") +
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

edr_differences <- to_plot

####### SD Comparison Plot ########################

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

sp <- c("LCTH", "LEPC", "HASP", "SPOW", "KIWA", "BITH")

to_plot <- dis_summary[which(dis_summary$Code %in% sp), ]

species_vars <- to_plot$variable

(species_prediction_plot <- bayesplot::mcmc_intervals(exp(dis_model$draws(species_vars)) * 100) + 
    xlab("Predicted EDR") + 
    ylab("Species") +
    scale_y_discrete(labels = (to_plot$Code)) +
    coord_flip())

####### Plotting Other Parameters #################

(covariates_factors <- bayesplot::mcmc_areas(dis_model$draws(c("intercept", 
                                                       "mu_mig_strat[1]",
                                                       "mu_mig_strat[2]",
                                                       "mu_habitat[1]",
                                                       "mu_habitat[2]")),
                                         prob = 0.95))

(covariates_slopes <- bayesplot::mcmc_areas(dis_model$draws(c("beta_mass", 
                                                               "beta_pitch")),
                                             prob = 0.95))

(sigma_plot <- bayesplot::mcmc_areas(dis_model$draws(c("sigma")),
                                     prob = 0.95))

params_summary <- dis_model$summary(c("intercept", "mu_mig_strat",
                                      "mu_habitat", "beta_mass", "beta_pitch",
                                      "sigma"))

####### Output ####################################

write.table(dis_summary, file = "data/generated/tau.csv", sep = ",", row.names = FALSE)
write.table(edr_differences, file = "data/generated/tau_differences.csv", sep = ",", row.names = FALSE)
write.table(params_summary, file = "data/generated/distance_params.csv", sep = ",", row.names = FALSE)

png("output/plots/distance_1vs1.png",
    width = 6, height = 3, res = 600, units = "in")
ggarrange(modelled_difference_plot, single_vs_multi_plot,
                    ncol = 2,
                    labels = c("A", "B"))
dev.off()

png("output/plots/distance_sd_plot.png",
    width = 6, height = 4, res = 600, units = "in")
ggarrange(sd_intercept_plot, sd_slope_plot, ncol = 2, labels = c("A", "B"))
dev.off()

png("output/plots/distance_predictions_plot.png",
    width = 6, height = 4, res = 600, units = "in")
print(species_prediction_plot)
dev.off()

png("output/plots/distance_other_params.png",
    width = 6, height = 6, res = 600, units = "in")
ggarrange(covariates_factors, covariates_slopes, sigma_plot, nrow = 3, labels = c("A", "B", "C"))
dev.off()
