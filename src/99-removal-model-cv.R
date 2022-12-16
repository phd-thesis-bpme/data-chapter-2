####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 99-removal-model-cv.R
# Created April 2022
# Last Updated December 2022

####### Import Libraries and External Files #######

library(cmdstanr)

####### Load Data #################################

load("data/generated/removal_stan_data_cv.rda")

####### Set Constants #############################

models <- c("brownian", "pagel")
removal_stan_data_cv$grainsize <- 1
removal_stan_data_cv$lambda <- 0.76

# Stan settings
n_iter <- 2000
n_warmup <- 1000
n_chains <- 4
refresh <- 100
threads_per_chain <- 4

####### Cross Validation ##########################

sp_list <- removal_stan_data_cv$sp_list
removal_stan_data_cv$sp_list <- NULL

for (m in models)
{
  dir.create(paste0("output/cv/removal/", m), recursive = TRUE)
  model_file <- cmdstan_model(stan_file = paste0("models/removal_", m, ".stan"),
                         cpp_options = list(stan_threads = TRUE))
  
  # first, run the full model with no species removed
  removal_stan_fit_full <- model_file$sample(
    data = removal_stan_data_cv,
    iter_warmup = n_warmup,
    iter_sampling = n_iter,
    chains = n_chains,
    parallel_chains = n_chains,
    refresh = refresh,
    threads_per_chain = threads_per_chain
  )
  removal_stan_fit_full$save_object(file = paste0("output/cv/removal/", m,"/1full.RDS"))
  rm(removal_stan_fit_fill)
  
  sp_list_unique <- unique(sp_list)[,1]
  
  for (sp in sp_list_unique)
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
      iter_warmup = n_warmup,
      iter_sampling = n_iter,
      chains = n_chains,
      parallel_chains = n_chains,
      refresh = refresh,
      threads_per_chain = threads_per_chain
    )
    removal_stan_fit_sp$save_object(file = paste0("output/cv/removal/",
                                                  m,
                                                  "/",
                                                  sp,
                                                  ".RDS"))
    rm(removal_stan_fit_sp)
  }
}
