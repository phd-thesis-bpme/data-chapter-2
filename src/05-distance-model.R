####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 05-distance-model.R
# Created October 2022
# Last Updated March 2023

####### Import Libraries and External Files #######

library(cmdstanr)
library(MASS)
library(ggplot2)
library(ggpubr)
theme_set(theme_pubclean())

####### Load Data #################################

load("data/generated/distance_stan_data_pred.rda")

####### Set Constants #############################

distance_stan_data_pred$grainsize <- 1
distance_stan_data_pred$lambda <- 0

# Scale the maximum distances to units of KM for computational ease
distance_stan_data_pred$max_dist <- distance_stan_data_pred$max_dist / 1000

# Stan settings
n_iter <- 2000
n_warmup <- 1000
n_chains <- 4
refresh <- 10
threads_per_chain <- 3

# Prior predictive check settings
n_sims <- 100

####### Prior Predictive Check ####################

phylo_corr_pl <- distance_stan_data_pred$phylo_corr * distance_stan_data_pred$lambda
for (i in 1:dim(phylo_corr_pl)[1])
{
  phylo_corr_pl[i,i] <- 1
}

# Sigma
sigma <- rexp(n = n_sims, rate = 5)
pdf(file = "output/prior_predictive_check/distance/sigma.pdf")
print(ggplot(data = data.frame(sigma), aes(x = sigma)) +
        geom_histogram(bins = 20) +
        NULL)
dev.off()

# mu mig strat
mu_mig_strat <- matrix(data = NA,
                       ncol = distance_stan_data_pred$n_mig_strat,
                       nrow = n_sims)
for (i in 1:distance_stan_data_pred$n_mig_strat)
{
  mu_mig_strat[,i] <- rnorm(n_sims, mean = 0, sd = 0.01)
}
pdf(file = "output/prior_predictive_check/distance/mu_mig_strat.pdf")
to_plot <- data.frame(Value = c(mu_mig_strat[,1],
                                mu_mig_strat[,2]),
                      Mig_Strat = c(rep("Resident", 100),
                                    rep("Migrant", 100)))
for (i in unique(to_plot$Mig_Strat))
{
  print(ggplot(data = to_plot[which(to_plot$Mig_Strat == i), ],
               aes(x = Value)) +
          geom_histogram(bins = 20) +
          xlab(i) +
          xlim(floor(min(to_plot$Value) - 1), ceiling(max(to_plot$Value) + 1)) +
          NULL)
}
dev.off()

# mu habitat
mu_habitat <- matrix(data = NA,
                       ncol = distance_stan_data_pred$n_habitat,
                       nrow = n_sims)
for (i in 1:distance_stan_data_pred$n_habitat)
{
  mu_habitat[,i] <- rnorm(n_sims, mean = 0, sd = 0.01)
}
pdf(file = "output/prior_predictive_check/distance/mu_habitat.pdf")
to_plot <- data.frame(Value = c(mu_habitat[,1],
                                mu_habitat[,2]),
                      Habitat = c(rep("Open", 100),
                                    rep("Closed", 100)))
for (i in unique(to_plot$Habitat))
{
  print(ggplot(data = to_plot[which(to_plot$Habitat == i), ],
               aes(x = Value)) +
          geom_histogram(bins = 20) +
          xlab(i) +
          xlim(floor(min(to_plot$Value) - 1), ceiling(max(to_plot$Value) + 1)) +
          NULL)
}
dev.off()

# beta mass
beta_mass <- rnorm(n = n_sims, mean = 0.01, sd = 0.01)
pdf(file = "output/prior_predictive_check/distance/beta_mass.pdf")
mass_plot <- ggplot() 
predictors <- seq(min(distance_stan_data_pred$mass),
                  max(distance_stan_data_pred$mass),
                  by = 0.01)
for (i in 1:n_sims)
{
  line_to_plot <- data.frame(Log_Mass = predictors,
                             Prediction = predictors * beta_mass[i])
  mass_plot <- mass_plot + 
    geom_line(data = line_to_plot, aes(x = Log_Mass, y = Prediction))
}
print(mass_plot)
dev.off()

# beta pitch
beta_pitch <- rnorm(n = n_sims, mean = -0.01, sd = 0.01)
pdf(file = "output/prior_predictive_check/distance/beta_pitch.pdf")
pitch_plot <- ggplot() 
predictors <- seq(min(distance_stan_data_pred$pitch),
                  max(distance_stan_data_pred$pitch),
                  by = 0.01)
for (i in 1:n_sims)
{
  line_to_plot <- data.frame(Log_Pitch = predictors,
                             Prediction = predictors * beta_pitch[i])
  pitch_plot <- pitch_plot + 
    geom_line(data = line_to_plot, aes(x = Log_Pitch, y = Prediction))
}
print(pitch_plot)
dev.off()

mu <- matrix(data = NA, nrow = n_sims, ncol = distance_stan_data_pred$n_species)
pdf(file = "output/prior_predictive_check/distance/mu.pdf")
for (s in 1:distance_stan_data_pred$n_species)
{
  mu[,s] <- (mu_mig_strat[, distance_stan_data_pred$mig_strat[s]]) +
                         (mu_habitat[, distance_stan_data_pred$habitat[s]]) +
                         (beta_mass * distance_stan_data_pred$mass[s]) +
                         (beta_pitch * distance_stan_data_pred$pitch[s])
  to_plot <- data.frame(Value = mu[,s])
  print(ggplot(data = to_plot, aes(x = Value)) +
          geom_histogram(bins = 20) +
          xlab(dimnames(phylo_corr_pl)[[1]][s]) +
          xlim(floor(min(mu) - 1), ceiling(max(mu) + 1)) +
          NULL)
}

dev.off()

log_tau <- matrix(data = NA, nrow = n_sims, ncol = distance_stan_data_pred$n_species)
for (i in 1:n_sims)
{
  log_tau[i,] <- MASS::mvrnorm(n = 1, mu = mu[i,],
                               Sigma = phylo_corr_pl * sigma[i])
}

pdf(file = "output/prior_predictive_check/distance/log_tau.pdf")
for (s in 1:distance_stan_data_pred$n_species)
{
  to_plot <- data.frame(Value = (log_tau[,s]))
  print(ggplot(data = to_plot, aes(x = Value)) +
          geom_histogram(bins = 20) +
          xlab(dimnames(phylo_corr_pl)[[1]][s]) +
          #xlim(floor(min(log_tau) - 1), ceiling(max(log_tau) + 1)) +
          NULL)
}
dev.off()

pdf(file = "output/prior_predictive_check/distance/tau.pdf")
tau <- exp(log_tau) * 1000
for (s in 1:distance_stan_data_pred$n_species)
{
  to_plot <- data.frame(Value = (tau[,s]))
  print(ggplot(data = to_plot, aes(x = Value)) +
          geom_histogram(bins = 20) +
          xlab(dimnames(phylo_corr_pl)[[1]][s]) +
          #xlim(0, 500) +
          NULL)
}
dev.off()

####### Run Model #################################

model_file <- cmdstan_model(stan_file = "models/distance.stan",
                            cpp_options = list(stan_threads = TRUE))

stan_run <- model_file$sample(
  data = distance_stan_data_pred,
  iter_warmup = n_warmup,
  iter_sampling = n_iter,
  chains = n_chains,
  parallel_chains = n_chains,
  refresh = refresh,
  threads_per_chain = threads_per_chain
)
stan_run$save_object(file = paste0("output/model_runs/distance_predictions.RDS"))

