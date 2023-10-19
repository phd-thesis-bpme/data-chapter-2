####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 4-distance-ms-vs-ss.R
# Created October 2023
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

####### Data Wrangling ############################

#' First we must drop all species which are just not going to be involved
#' in the cross validation at all. This will be the species for which we are
#' generated predictions
pred_drops <- c("LCTH", "LEPC", "HASP", "SPOW", "KIWA", "BITH")
dis_data <- subset_distance_data(distance_stan_data = distance_stan_data,
                                 sps = setdiff(unique(distance_stan_data$sp_list),
                                               pred_drops))
rm(distance_stan_data)

ms_model <- cmdstan_model(stan_file = "models/distance_cp.stan",
                          cpp_options = list(stan_threads = TRUE))

ss_model <- cmdstan_model(stan_file = "models/distance_single.stan",
                          cpp_options = list(stan_threads = TRUE))

####### Run Full Models ###########################

inits <- generate_distance_inits(n_chains = n_chains,
                                 sp_list = setdiff(as.vector(dis_data$sp_all), pred_drops),
                                 napops_skip = NULL,
                                 param = "cp")

# Multi species model
stan_run_ms <- ms_model$sample(
  data = dis_data,
  iter_warmup = n_warmup,
  iter_sampling = n_iter,
  chains = n_chains,
  parallel_chains = n_chains,
  refresh = refresh,
  threads_per_chain = threads_per_chain,
  init = inits
)
stan_run_ms$save_object(file = paste0("output/model_runs/distance_ms.RDS"))

# Single species model
dis_data2 <- dis_data
dis_data2$n_mig_strat <- NULL
dis_data2$mig_strat <- NULL
dis_data2$n_habitat <- NULL
dis_data2$habitat <- NULL
dis_data2$mass <- NULL
dis_data2$pitch <- NULL

stan_run_ss <- ss_model$sample(
  data = dis_data2,
  iter_warmup = n_warmup,
  iter_sampling = n_iter,
  chains = n_chains,
  parallel_chains = n_chains,
  refresh = refresh,
  threads_per_chain = threads_per_chain,
  init = inits
)
stan_run_ss$save_object(file = paste0("output/model_runs/distance_ss.RDS"))
