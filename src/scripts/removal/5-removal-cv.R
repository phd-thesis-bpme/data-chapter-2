####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 5-removal-cv.R
# Created October 2023
# Last Updated March 2024

####### Import Libraries and External Files #######

library(cmdstanr)

####### Load Data #################################

load("data/generated/removal_stan_data_cv.rda")
cv_folds <- read.csv("data/generated/removal_cv_folds.csv")

####### Set Constants #############################

# Stan settings
n_iter <- 500
n_warmup <- 1000
n_chains <- 4
refresh <- 10
threads_per_chain <- 3

####### Data Wrangling ############################

ms_model <- cmdstan_model(stan_file = "models/removal.stan",
                          cpp_options = list(stan_threads = TRUE))

ss_model <- cmdstan_model(stan_file = "models/removal_single.stan",
                          cpp_options = list(stan_threads = TRUE))

####### Run Cross Validation ######################

for (i in 1:max(cv_folds$cv_fold))
{
  indices_to_drop <- which(cv_folds$cv_fold == i)
  
  stan_cv_data <- list(n_samples = removal_stan_data_cv$n_samples - length(indices_to_drop),
                       n_species = removal_stan_data_cv$n_species,
                       max_intervals = removal_stan_data_cv$max_intervals,
                       species = removal_stan_data_cv$species[-indices_to_drop],
                       abund_per_band = removal_stan_data_cv$abund_per_band[-indices_to_drop, ],
                       bands_per_sample = removal_stan_data_cv$bands_per_sample[-indices_to_drop],
                       max_time = removal_stan_data_cv$max_time[-indices_to_drop, ],
                       phylo_corr = removal_stan_data_cv$phylo_corr,
                       n_mig_strat = removal_stan_data_cv$n_mig_strat,
                       mig_strat = removal_stan_data_cv$mig_strat,
                       sp_all = removal_stan_data_cv$sp_all)
  
  # Check to see if we can drop any training columns (i.e. if bands per sample is smaller)
  if (max(stan_cv_data$bands_per_sample) < stan_cv_data$max_intervals)
  {
    max_intervals <- max(stan_cv_data$bands_per_sample)
    
    stan_cv_data$max_intervals <- max_intervals
    stan_cv_data$abund_per_band <- stan_cv_data$abund_per_band[, 1:max_intervals]
    stan_cv_data$max_time <- stan_cv_data$max_time[, 1:max_intervals]
  }
  
  stan_cv_data$grainsize <- 1
  stan_cv_data$lambda <- 0.79
  
  stan_cv_data$sp_all <- NULL
  
  stan_run_ms <- ms_model$sample(
    data = stan_cv_data,
    iter_warmup = n_warmup,
    iter_sampling = n_iter,
    chains = n_chains,
    parallel_chains = n_chains,
    refresh = refresh,
    threads_per_chain = threads_per_chain
  )
  stan_run_ms$save_object(file = paste0("output/model_runs/cv_removal/ms-vs-ss/ms_fold_",
                                        i,
                                        ".RDS"))
  
  stan_cv_data$n_mig_strat <- NULL
  stan_cv_data$mig_strat <- NULL
  stan_cv_data$phylo_corr <- NULL
  stan_cv_data$lambda <- NULL
  
  stan_run_ss <- ss_model$sample(
    data = stan_cv_data,
    iter_warmup = n_warmup,
    iter_sampling = n_iter,
    chains = n_chains,
    parallel_chains = n_chains,
    refresh = refresh,
    threads_per_chain = threads_per_chain
  )
  stan_run_ss$save_object(file = paste0("output/model_runs/cv_removal/ms-vs-ss/ss_fold_",
                                        i,
                                        ".RDS"))
}
