####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 02-removal-model-cv.R
# Created April 2022
# Last Updated October 2022

####### Import Libraries and External Files #######

library(cmdstanr)

####### Load Data #################################

load("data/generated/removal_stan_data_cv.rda")

####### Cross Validation ##########################

sp_list <- removal_stan_data_cv$sp_list
removal_stan_data_cv$sp_list <- NULL
models <- c("brownian", "ou")

for (m in models)
{
  dir.create(paste0("output/cv/removal/", m), recursive = TRUE)
  model_file <- cmdstan_model(stan_file = paste0("models/removal_", m, ".stan"),
                         cpp_options = list(stan_threads = TRUE))
  removal_stan_data_cv$grainsize <- 1
  
  # first, run the full model with no species removed
  removal_stan_fit_full <- model_file$sample
  (
    data = removal_stan_data,
    iter_warmup = 10,
    iter_sampling = 50,
    chains = 4,
    parallel_chains = 4,
    refresh = 10,
    threads_per_chain = 4
  )
  removal_stan_fit_full$save_object(file = paste0("output/cv/removal/", m,"/1full.RDS"))
  
}

####### Output ####################################


#save(removal_stan_fit, file = "data/generated/removal_model.rda")

