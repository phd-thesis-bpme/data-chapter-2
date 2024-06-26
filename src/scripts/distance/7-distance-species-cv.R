####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 7-distance-species-cv.R
# Created March 2024
# Last Updated May 2024

####### Import Libraries and External Files #######

library(cmdstanr)

source("src/functions/generate-distance-inits.R")

####### Load Data #################################

load("data/generated/distance_stan_data_cv.rda")
cv_folds <- read.csv("data/generated/distance_cv_folds_species.csv")

####### Set Constants #############################

# Stan settings
n_iter <- 2000
n_warmup <- 1000
n_chains <- 4
refresh <- 10
threads_per_chain <- 7

####### Data Wrangling ############################

model <- cmdstan_model(stan_file = "models/distance_cp.stan",
                          cpp_options = list(stan_threads = TRUE))

####### Run Cross Validation ######################

for (i in 1:max(cv_folds$cv_fold))
{
  indices_to_drop <- which(cv_folds$cv_fold == i)
  
  stan_cv_data <- list(n_samples = distance_stan_data_cv$n_samples - length(indices_to_drop),
                       n_species = distance_stan_data_cv$n_species,
                       max_intervals = distance_stan_data_cv$max_intervals,
                       species = distance_stan_data_cv$species[-indices_to_drop],
                       abund_per_band = distance_stan_data_cv$abund_per_band[-indices_to_drop, ],
                       bands_per_sample = distance_stan_data_cv$bands_per_sample[-indices_to_drop],
                       max_dist = distance_stan_data_cv$max_dist[-indices_to_drop, ],
                       n_mig_strat = distance_stan_data_cv$n_mig_strat,
                       mig_strat = distance_stan_data_cv$mig_strat,
                       n_habitat = distance_stan_data_cv$n_habitat,
                       habitat = distance_stan_data_cv$habitat,
                       mass = distance_stan_data_cv$mass,
                       pitch = distance_stan_data_cv$pitch,
                       sp_all = distance_stan_data_cv$sp_all)
  
  # Check to see if we can drop any training columns (i.e. if bands per sample is smaller)
  if (max(stan_cv_data$bands_per_sample) < stan_cv_data$max_intervals)
  {
    max_intervals <- max(stan_cv_data$bands_per_sample)
    
    stan_cv_data$max_intervals <- max_intervals
    stan_cv_data$abund_per_band <- stan_cv_data$abund_per_band[, 1:max_intervals]
    stan_cv_data$max_dist <- stan_cv_data$max_dist[, 1:max_intervals]
  }
  
  stan_cv_data$grainsize <- 1
  stan_cv_data$max_dist <- stan_cv_data$max_dist / 100
  
  inits <- generate_distance_inits(n_chains = n_chains,
                                   sp_list = stan_cv_data$sp_all,
                                   napops_skip = NULL,
                                   param = "cp")
  stan_cv_data$sp_all <- NULL
  
  # Strip away center and scale attributes from mass and pitch for analysis
  stan_cv_data$mass <- stan_cv_data$mass[,1]
  stan_cv_data$pitch <- stan_cv_data$pitch[,1]
  
  stan_run <- model$sample(
    data = stan_cv_data,
    iter_warmup = n_warmup,
    iter_sampling = n_iter,
    chains = n_chains,
    parallel_chains = n_chains,
    refresh = refresh,
    threads_per_chain = threads_per_chain,
    init = inits
  )
  stan_run$save_object(file = paste0("output/model_runs/cv_distance/species/fold_",
                                        i,
                                        ".RDS"))
}
