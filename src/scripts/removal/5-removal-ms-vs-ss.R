####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 5-removal-ms-vs-ss.R
# Created October 2023
# Last Updated March 2024

####### Import Libraries and External Files #######

library(cmdstanr)
library(dplyr)

####### Load Data #################################

load("data/generated/removal_stan_data_cv.rda")

####### Set Constants #############################

# Stan settings
n_iter <- 2000
n_warmup <- 1000
n_chains <- 4
refresh <- 10
threads_per_chain <- 4

####### Data Wrangling ############################

ms_model <- cmdstan_model(stan_file = "models/removal.stan",
                          cpp_options = list(stan_threads = TRUE))

ss_model <- cmdstan_model(stan_file = "models/removal_single.stan",
                          cpp_options = list(stan_threads = TRUE))

####### Run Full Models ###########################

removal_stan_data_cv$sp_list <- NULL
removal_stan_data_cv$sp_all <- NULL
removal_stan_data_cv$grainsize <- 1
removal_stan_data_cv$lambda <- 0.79

# Multi species model
stan_run_ms <- ms_model$sample(
  data = removal_stan_data_cv,
  iter_warmup = n_warmup,
  iter_sampling = n_iter,
  chains = n_chains,
  parallel_chains = n_chains,
  refresh = refresh,
  threads_per_chain = threads_per_chain
)
stan_run_ms$save_object(file = paste0("output/model_runs/removal_ms.RDS"))

# Single species model
removal_stan_data_cv$phylo_corr <- NULL
removal_stan_data_cv$m_mig_strat <- NULL
removal_stan_data_cv$mig_strat <- NULL

stan_run_ss <- ss_model$sample(
  data = removal_stan_data_cv,
  iter_warmup = n_warmup,
  iter_sampling = n_iter,
  chains = n_chains,
  parallel_chains = n_chains,
  refresh = refresh,
  threads_per_chain = threads_per_chain
)
stan_run_ss$save_object(file = paste0("output/model_runs/removal_ss.RDS"))
