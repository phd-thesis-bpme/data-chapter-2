####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 6-distance-predictions.R
# Created October 2022
# Last Updated October 2023

####### Import Libraries and External Files #######

library(cmdstanr)
library(dplyr)

source("src/functions/generate-distance-inits.R")
source("src/functions/subset-distance-data.R")

####### Load Data #################################

load("data/generated/distance_stan_data.rda")

####### Set Constants #############################

# Stan settings
n_iter <- 2000
n_warmup <- 1000
n_chains <- 4
refresh <- 10
threads_per_chain <- 7

distance_stan_data$grainsize <- 1

inits <- generate_distance_inits(n_chains = n_chains,
                                 napops_skip = c("BITH", "HASP", "KIWA", "LCTH", "LEPC", "SPOW"),
                                 sp_list = distance_stan_data$sp_all,
                                 param = "cp",
                                 species_cp = distance_stan_data$species_cp,
                                 species_ncp = distance_stan_data$species_ncp)

distance_stan_data$sp_list <- NULL
distance_stan_data$sp_all <- NULL

# Scale the maximum distances for computational ease
distance_stan_data$max_dist <- distance_stan_data$max_dist / 100

# get rid of centre/scale attributes for modelling
distance_stan_data$pitch <- distance_stan_data$pitch[,1]
distance_stan_data$mass <- distance_stan_data$mass[,1]

####### Run Model #################################

model_file <- cmdstan_model(stan_file = "models/distance_cp.stan",
                            cpp_options = list(stan_threads = TRUE))

stan_run <- model_file$sample(
  data = distance_stan_data,
  iter_warmup = n_warmup,
  iter_sampling = n_iter,
  chains = n_chains,
  parallel_chains = n_chains,
  refresh = refresh,
  threads_per_chain = threads_per_chain,
  output_dir = "output/model_runs/stan_output/",
  init = inits
)

stan_run$save_object(file = paste0("output/model_runs/distance_predictions.RDS"))
