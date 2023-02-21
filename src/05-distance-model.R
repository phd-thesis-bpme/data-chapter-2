####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 05-distance-model.R
# Created October 2022
# Last Updated February 2023

####### Import Libraries and External Files #######

library(cmdstanr)

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

