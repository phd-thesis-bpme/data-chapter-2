####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 03-removal-model.R
# Created December 2022
# Last Updated December 2022

####### Import Libraries and External Files #######

library(cmdstanr)

####### Load Data #################################

load("data/generated/removal_stan_data_pred.rda")

####### Set Constants #############################

models <- c("pagel")
removal_stan_data_cv$grainsize <- 1
removal_stan_data_cv$lambda <- 0.79

# Stan settings
n_iter <- 2000
n_warmup <- 1000
n_chains <- 4
refresh <- 100
threads_per_chain <- 3

####### Run Model #################################

dir.create(paste0("output/models/removal/"), recursive = TRUE)
model_file <- cmdstan_model(stan_file = "models/removal_pagel.stan",
                            cpp_options = list(stan_threads = TRUE))

stan_run <- model_file$sample(
  data = removal_stan_data_cv,
  iter_warmup = n_warmup,
  iter_sampling = n_iter,
  chains = n_chains,
  parallel_chains = n_chains,
  refresh = refresh,
  threads_per_chain = threads_per_chain
)
stan_run$save_object(file = paste0("output/cv/removal/", m,"/1full.RDS"))

