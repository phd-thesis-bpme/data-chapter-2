####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 05-distance-prior-check.R
# Created March 2023
# Last Updated May 2023

####### Import Libraries and External Files #######

library(MASS)
library(ggplot2)
library(ggpubr)
theme_set(theme_pubclean())

####### Load Data #################################

load("data/generated/distance_stan_data.rda")

####### Set Constants #############################

distance_stan_data$lambda <- 0

# Scale the maximum distances to units of KM for computational ease
#distance_stan_data$max_dist <- distance_stan_data$max_dist / 1000

# Prior predictive check settings
n_sims <- 100

####### Prior Predictive Check ####################

phylo_corr_pl <- distance_stan_data$phylo_corr * distance_stan_data$lambda
for (i in 1:dim(phylo_corr_pl)[1])
{
  phylo_corr_pl[i,i] <- 1
}

# Sigma
sigma <- rexp(n = n_sims, rate = 5)
pdf(file = "output/prior_predictive_check/distance/01-sigma.pdf")
print(ggplot(data = data.frame(sigma), aes(x = sigma)) +
        geom_histogram(bins = 20) +
        NULL)
dev.off()

intercept <- rnorm(n = n_sims, mean = 0.05, sd = 0.1)
pdf(file = "output/prior_predictive_check/distance/02-intercept.pdf")
print(ggplot(data = data.frame(intercept), aes(x = intercept)) +
        geom_histogram(bins = 20) +
        NULL)
dev.off()

# mu mig strat
mu_mig_strat <- matrix(data = NA,
                       ncol = distance_stan_data$n_mig_strat,
                       nrow = n_sims)

mu_mig_strat[,1] <- 0
mu_mig_strat[,2] <- rnorm(n_sims, sd = 0.05)
# for (i in 1:distance_stan_data$n_mig_strat)
# {
#   mu_mig_strat[,i] <- rnorm(n_sims)
# }
pdf(file = "output/prior_predictive_check/distance/03-mu_mig_strat.pdf")
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
          #xlim(floor(min(to_plot$Value) - 1), ceiling(max(to_plot$Value) + 1)) +
          NULL)
}
dev.off()

# mu habitat
mu_habitat <- matrix(data = NA,
                     ncol = distance_stan_data$n_habitat,
                     nrow = n_sims)

mu_habitat[,1] <- 0
mu_habitat[,2] <- rnorm(n_sims, sd = 0.05)
# for (i in 1:distance_stan_data$n_habitat)
# {
#   mu_habitat[,i] <- rnorm(n_sims)
# }
pdf(file = "output/prior_predictive_check/distance/04-mu_habitat.pdf")
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
          #xlim(floor(min(to_plot$Value) - 1), ceiling(max(to_plot$Value) + 1)) +
          NULL)
}
dev.off()

# beta mass
beta_mass <- rnorm(n = n_sims, mean = 0.01, sd = 0.005)
pdf(file = "output/prior_predictive_check/distance/05-beta_mass.pdf")

mass_hist <- ggplot(data = data.frame(Mass_Slope = beta_mass), aes(x = Mass_Slope)) +
  geom_histogram(bins = 20) +
  NULL
print(mass_hist)

mass_plot <- ggplot() 
predictors <- seq(min(distance_stan_data$mass),
                  max(distance_stan_data$mass),
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
beta_pitch <- rnorm(n = n_sims, mean = -0.01, sd = 0.005)
pdf(file = "output/prior_predictive_check/distance/06-beta_pitch.pdf")

pitch_hist <- ggplot(data = data.frame(Pitch_Slope = beta_pitch), aes(x = Pitch_Slope)) +
  geom_histogram(bins = 20) +
  NULL
print(pitch_hist)

pitch_plot <- ggplot() 
predictors <- seq(min(distance_stan_data$pitch),
                  max(distance_stan_data$pitch),
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

mu <- matrix(data = NA, nrow = n_sims, ncol = distance_stan_data$n_species)
pdf(file = "output/prior_predictive_check/distance/07-mu.pdf")
for (s in 1:distance_stan_data$n_species)
{
  mu[,s] <- intercept + (mu_mig_strat[, distance_stan_data$mig_strat[s]]) +
    (mu_habitat[, distance_stan_data$habitat[s]]) +
    (beta_mass * distance_stan_data$mass[s]) +
    (beta_pitch * distance_stan_data$pitch[s])
  to_plot <- data.frame(Value = mu[,s])
  print(ggplot(data = to_plot, aes(x = Value)) +
          geom_histogram(bins = 20) +
          xlab(dimnames(phylo_corr_pl)[[1]][s]) +
          xlim(floor(min(mu) - 1), ceiling(max(mu) + 1)) +
          NULL)
}

dev.off()

log_tau <- matrix(data = NA, nrow = n_sims, ncol = distance_stan_data$n_species)
for (i in 1:n_sims)
{
  log_tau[i,] <- MASS::mvrnorm(n = 1, mu = mu[i,],
                               Sigma = phylo_corr_pl * sigma[i])
}

pdf(file = "output/prior_predictive_check/distance/08-log_tau.pdf")
for (s in 1:distance_stan_data$n_species)
{
  to_plot <- data.frame(Value = (log_tau[,s]))
  print(ggplot(data = to_plot, aes(x = Value)) +
          geom_histogram(bins = 20) +
          xlab(dimnames(phylo_corr_pl)[[1]][s]) +
          #xlim(floor(min(log_tau) - 1), ceiling(max(log_tau) + 1)) +
          NULL)
}
dev.off()

pdf(file = "output/prior_predictive_check/distance/09-tau.pdf")
tau <- exp(log_tau) * 100
for (s in 1:distance_stan_data$n_species)
{
  to_plot <- data.frame(Value = (tau[,s]))
  print(ggplot(data = to_plot, aes(x = Value)) +
          geom_histogram(bins = 20) +
          xlab(dimnames(phylo_corr_pl)[[1]][s]) +
          #xlim(0, 500) +
          NULL)
}
dev.off()
