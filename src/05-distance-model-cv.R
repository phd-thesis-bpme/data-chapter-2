####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 05-distance-model-cv.R
# Created October 2022
# Last Updated October 2022

####### Import Libraries and External Files #######

library(cmdstanr)

####### Load Data #################################

load("data/generated/distance_stan_data_cv.rda")

####### Set Constants #############################

models <- c("brownian", "pagel")
distance_stan_data_cv$grainsize <- 1
distance_stan_data_cv$lambda <- 0.76

# Stan settings
n_iter <- 2000
n_warmup <- 1000
n_chains <- 4
refresh <- 100
threads_per_chain <- 4

####### Cross Validation ##########################

sp_list <- distance_stan_data_cv$sp_list
distance_stan_data_cv$sp_list <- NULL

for (m in models)
{
  dir.create(paste0("output/cv/distance/", m), recursive = TRUE)
  model_file <- cmdstan_model(stan_file = paste0("models/distance_", m, ".stan"),
                         cpp_options = list(stan_threads = TRUE))
  
  # first, run the full model with no species removed
  distance_stan_fit_full <- model_file$sample(
    data = distance_stan_data_cv,
    iter_warmup = n_warmup,
    iter_sampling = n_iter,
    chains = n_chains,
    parallel_chains = n_chains,
    refresh = refresh,
    threads_per_chain = threads_per_chain
  )
  distance_stan_fit_full$save_object(file = paste0("output/cv/distance/", m,"/1full.RDS"))
  rm(distance_stan_fit_fill)
  
  sp_list_unique <- unique(sp_list)[,1]
  
  for (sp in sp_list_unique)
  {
    #' Figure out which indices need to be removed, and create a new data list
    #' with those particular indices removed from each matrix or vector.
    #' However, we want to keep the number of species the same, as well
    #' as the correlation matrix, because we need to try to predict for the
    #' species which now does not have data.

    indices_to_remove <- which(sp_list == sp)
    distance_stan_data_cv_sp <- distance_stan_data_cv
    distance_stan_data_cv_sp$species <- 
      distance_stan_data_cv_sp$species[-c(indices_to_remove)]
    distance_stan_data_cv_sp$abund_per_band <- 
      distance_stan_data_cv_sp$abund_per_band[-c(indices_to_remove),]
    distance_stan_data_cv_sp$bands_per_sample <- 
      distance_stan_data_cv_sp$bands_per_sample[-c(indices_to_remove)]
    distance_stan_data_cv_sp$max_time <- 
      distance_stan_data_cv_sp$max_time[-c(indices_to_remove),]
    distance_stan_data_cv_sp$n_samples <- 
      length(distance_stan_data_cv_sp$species)
    
    distance_stan_fit_sp <- model_file$sample(
      data = distance_stan_data_cv_sp,
      iter_warmup = n_warmup,
      iter_sampling = n_iter,
      chains = n_chains,
      parallel_chains = n_chains,
      refresh = refresh,
      threads_per_chain = threads_per_chain
    )
    distance_stan_fit_sp$save_object(file = paste0("output/cv/distance/",
                                                  m,
                                                  "/",
                                                  sp,
                                                  ".RDS"))
    rm(distance_stan_fit_sp)
  }
}
