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
models <- c("ou")
removal_stan_data_cv$grainsize <- 1

for (m in models)
{
  dir.create(paste0("output/cv/removal/", m), recursive = TRUE)
  model_file <- cmdstan_model(stan_file = paste0("models/removal_", m, ".stan"),
                         cpp_options = list(stan_threads = TRUE))

  
  # first, run the full model with no species removed
  removal_stan_fit_full <- model_file$sample(
    data = removal_stan_data_cv,
    iter_warmup = 10,
    iter_sampling = 50,
    chains = 4,
    parallel_chains = 4,
    refresh = 10,
    threads_per_chain = 4
  )
  removal_stan_fit_full$save_object(file = paste0("output/cv/removal/", m,"/1full.RDS"))
  
  sp_list_unique <- unique(sp_list)[,1]
  
  for (sp in sp_list_unique[1:3])
  {
    #' Figure out which indices need to be removed, and create a new data list
    #' with those particular indices removed from each matrix or vector.
    #' However, we want to keep the number of species the same, as well
    #' as the correlation matrix, because we need to try to predict for the
    #' species which now does not have data.

    indices_to_remove <- which(sp_list == sp)
    removal_stan_data_cv_sp <- removal_stan_data_cv
    removal_stan_data_cv_sp$species <- 
      removal_stan_data_cv_sp$species[-c(indices_to_remove)]
    removal_stan_data_cv_sp$abund_per_band <- 
      removal_stan_data_cv_sp$abund_per_band[-c(indices_to_remove),]
    removal_stan_data_cv_sp$bands_per_sample <- 
      removal_stan_data_cv_sp$bands_per_sample[-c(indices_to_remove)]
    removal_stan_data_cv_sp$max_time <- 
      removal_stan_data_cv_sp$max_time[-c(indices_to_remove),]
    removal_stan_data_cv_sp$n_samples <- 
      length(removal_stan_data_cv_sp$species)
    
    removal_stan_fit_sp <- model_file$sample(
      data = removal_stan_data_cv_sp,
      iter_warmup = 10,
      iter_sampling = 50,
      chains = 4,
      parallel_chains = 4,
      refresh = 10,
      threads_per_chain = 4
    )
    removal_stan_fit_sp$save_object(file = paste0("output/cv/removal/",
                                                  m,
                                                  "/",
                                                  sp,
                                                  ".RDS"))
  }
}
