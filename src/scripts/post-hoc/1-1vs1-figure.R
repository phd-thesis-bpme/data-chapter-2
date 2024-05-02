####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# posthoc/1-1vs1-figure.R
# Created December 2023
# Last Updated May 2024

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

rem_ms <- readRDS("output/model_runs/removal_ms.RDS")
rem_ss <- readRDS("output/model_runs/removal_ss.RDS")
load("data/generated/removal_stan_data_cv.rda")

dis_ms <- readRDS("output/model_runs/distance_ms.RDS")
dis_ss <- readRDS("output/model_runs/distance_ss.RDS")
load("data/generated/distance_stan_data_cv.rda")

####### Removal Model #############################

# Extract log_phi summary statistics from full Stan model runs and add species Codes
rem_ms_summary <- rem_ms$summary("log_phi")
rem_ms_summary$Code <- removal_stan_data_cv$sp_all

rem_ss_summary <- rem_ss$summary("log_phi")
rem_ss_summary$Code <- removal_stan_data_cv$sp_all

# Get data sample size for all species and add to summary
species_n <- data.frame(table(removal_stan_data_cv$species))
names(species_n) <- c("index", "N")

rem_ms_summary$index <- seq(1, nrow(rem_ms_summary))
rem_ms_summary$N <- 0
rem_ss_summary$index <- seq(1, nrow(rem_ss_summary))
rem_ss_summary$N <- 0
for (i in 1:nrow(species_n))
{
  rem_ms_summary[which(rem_ms_summary$index == species_n$index[i]), "N"] <-
    species_n$N[i]
  
  rem_ss_summary[which(rem_ss_summary$index == species_n$index[i]), "N"] <-
    species_n$N[i]
}

####### 1-to-1 Plot ###############################

to_plot <- data.frame(Species = rem_ms_summary$Code,
                      Single = rem_ss_summary$mean,
                      Multi = rem_ms_summary$mean,
                      N = rem_ms_summary$N)
to_plot$Label <- ""
to_plot$Multi <- exp(to_plot$Multi)
to_plot$Single <- exp(to_plot$Single)
to_plot$diff <- to_plot$Multi - to_plot$Single

# Check for a 10% relative difference in cue rates
for (i in 1:nrow(to_plot))
{
  if ((abs(to_plot$diff[i]) / to_plot$Single[i]) > 0.10)
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

cr_differences <- to_plot
cr_diff_model <- diff_model_run$summary()

####### Distance Model ############################

# Extract log_tau summary statistics from full Stan model runs and add species Codes
dis_ms_summary <- dis_ms$summary("log_tau")
dis_ms_summary$Code <- distance_stan_data_cv$sp_all

dis_ss_summary <- dis_ss$summary("log_tau")
dis_ss_summary$Code <- distance_stan_data_cv$sp_all

# Get data sample size for all species and add to summary
species_n <- data.frame(table(distance_stan_data_cv$species))
names(species_n) <- c("index", "N")

dis_ms_summary$index <- seq(1, nrow(dis_ms_summary))
dis_ms_summary$N <- 0
dis_ss_summary$index <- seq(1, nrow(dis_ss_summary))
dis_ss_summary$N <- 0
for (i in 1:nrow(species_n))
{
  dis_ms_summary[which(dis_ms_summary$index == species_n$index[i]), "N"] <-
    species_n$N[i]
  
  dis_ss_summary[which(dis_ss_summary$index == species_n$index[i]), "N"] <-
    species_n$N[i]
}

to_plot <- data.frame(Species = dis_ms_summary$Code,
                      Single = dis_ss_summary$mean,
                      Multi = dis_ms_summary$mean,
                      N = dis_ms_summary$N)
to_plot$Label <- ""
to_plot$Multi <- exp(to_plot$Multi) * 100
to_plot$Single <- exp(to_plot$Single) * 100
to_plot$diff <- to_plot$Multi - to_plot$Single

# Check for a 20% relative difference in cue rates
for (i in 1:nrow(to_plot))
{
  if ((abs(to_plot$diff[i]) / to_plot$Single[i]) > 0.10)
  {
    to_plot$Label[i] <- to_plot$Species[i]
  }
}

(single_vs_multi_plot <- ggplot(data = to_plot, mapping = aes(x = Single, y = Multi)) +
    geom_point() +
    geom_abline(slope = 1) +
    xlim(0,800) +
    ylim(0,800) +
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

edr_differences <- to_plot
edr_diff_model <- diff_model_run$summary()

####### Output ####################################

write.table(rem_ms_summary, file = "data/generated/phi_ms.csv", sep = ",", row.names = FALSE)
write.table(rem_ss_summary, file = "data/generated/phi_ss.csv", sep = ",", row.names = FALSE)
write.table(dis_ms_summary, file = "data/generated/tau_ms.csv", sep = ",", row.names = FALSE)
write.table(dis_ss_summary, file = "data/generated/tau_ss.csv", sep = ",", row.names = FALSE)

write.table(cr_differences, file = "data/generated/phi_differences.csv", sep = ",", row.names = FALSE)
write.table(cr_diff_model, file = "data/generated/phi_difference_model.csv", sep = ",", row.names = FALSE)
write.table(edr_differences, file = "data/generated/tau_differences.csv", sep = ",", row.names = FALSE)
write.table(edr_diff_model, file = "data/generated/tau_difference_model.csv", sep = ",", row.names = FALSE)


tiff("output/plots/1vs1.tiff",
     width = 6, height = 6, res = 600, units = "in")
ggarrange(removal_plot, distance_plot, nrow = 2, labels = c("A", "B"))
dev.off()

png("output/plots/1vs1.png",
     width = 6, height = 6, res = 600, units = "in")
ggarrange(removal_plot, distance_plot, nrow = 2, labels = c("A", "B"))
dev.off()
