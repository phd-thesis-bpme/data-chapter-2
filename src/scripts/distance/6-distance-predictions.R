####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 6-distance-predictions.R
# Created October 2022
# Last Updated April 2024

####### Import Libraries and External Files #######

library(cmdstanr)

source("src/functions/generate-distance-inits.R")

####### Load Data #################################

load("data/generated/distance_stan_data_pred.rda")

####### Set Constants #############################

# Stan settings
n_iter <- 2000
n_warmup <- 1000
n_chains <- 4
refresh <- 10
threads_per_chain <- 4

distance_stan_data_pred$grainsize <- 1

inits <- generate_distance_inits(n_chains = n_chains,
                                 napops_skip = c("BITH", "HASP", "KIWA", "LCTH", "LEPC", "SPOW"),
                                 sp_list = distance_stan_data_pred$sp_all,
                                 param = "cp")

distance_stan_data_pred$sp_list <- NULL
distance_stan_data_pred$sp_all <- NULL

# Scale the maximum distances for computational ease
distance_stan_data_pred$max_dist <- distance_stan_data_pred$max_dist / 100

# get rid of centre/scale attributes for modelling
distance_stan_data_pred$pitch <- distance_stan_data_pred$pitch[,1]
distance_stan_data_pred$mass <- distance_stan_data_pred$mass[,1]

####### Run Model #################################

model_file <- cmdstan_model(stan_file = "models/distance_cp.stan",
                            cpp_options = list(stan_threads = TRUE))

stan_run <- model_file$sample(
  data = distance_stan_data_pred,
  iter_warmup = n_warmup,
  iter_sampling = n_iter,
  chains = n_chains,
  parallel_chains = n_chains,
  refresh = refresh,
  threads_per_chain = threads_per_chain,
  init = inits
)

stan_run$save_object(file = paste0("output/model_runs/distance_predictions.RDS"))
