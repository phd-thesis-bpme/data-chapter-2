####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 07-removal-model.R
# Created December 2022
# Last Updated October 2023

####### Import Libraries and External Files #######

library(cmdstanr)

####### Load Data #################################

load("data/generated/removal_stan_data.rda")

####### Set Constants #############################

removal_stan_data$grainsize <- 1
removal_stan_data$lambda <- 0.79

# Stan settings
n_iter <- 2000
n_warmup <- 1000
n_chains <- 4
refresh <- 10
threads_per_chain <- 3

removal_stan_data$sp_list <- NULL
removal_stan_data$sp_all <- NULL

####### Run Model #################################

model_file <- cmdstan_model(stan_file = "models/removal.stan",
                            cpp_options = list(stan_threads = TRUE))

stan_run <- model_file$sample(
  data = removal_stan_data,
  iter_warmup = n_warmup,
  iter_sampling = n_iter,
  chains = n_chains,
  parallel_chains = n_chains,
  refresh = refresh,
  threads_per_chain = threads_per_chain
)
stan_run$save_object(file = paste0("output/model_runs/removal_predictions.RDS"))

