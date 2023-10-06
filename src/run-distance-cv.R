####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 07-run-distance-cv.R
# Created October 2023
# Last Updated October 2023

####### Import Libraries and External Files #######

library(cmdstanr)
library(dplyr)

source("src/functions/generate-distance-inits.R")
source("src/functions/subset-distance-data.R")

####### Load Data #################################

load("data/generated/distance_stan_data.rda")
cv_folds <- read.csv("data/generated/distance_cv_folds.csv")

####### Set Constants #############################

# Stan settings
n_iter <- 2000
n_warmup <- 1000
n_chains <- 4
refresh <- 10
threads_per_chain <- 7

####### Run Cross Validation ######################

#' First we must drop all species which are just not going to be involved
#' in the cross validation at all. This will be the species for which we are
#' generated predictions
pred_drops <- c("LCTH", "LEPC", "HASP", "SPOW", "KIWA", "BITH")
dis_data <- subset_distance_data(distance_stan_data = distance_stan_data,
                                 sps = setdiff(unique(distance_stan_data$sp_list),
                                               pred_drops))
rm(distance_stan_data)

model_file <- cmdstan_model(stan_file = "models/distance_cp.stan",
                            cpp_options = list(stan_threads = TRUE))

for (i in 1:max(cv_folds$cv_fold))
{
  drops <- cv_folds[which(cv_folds$cv_fold == i), "Species"]
  indices_to_drop <- which(dis_data$sp_list %in% drops)
  
}




distance_stan_data$grainsize <- 1
distance_stan_data$lambda <- 0

distance_stan_data$sp_list <- NULL
phylo_corr <- distance_stan_data$phylo_corr
distance_stan_data$phylo_corr <- NULL

# Scale the maximum distances for computational ease
distance_stan_data$max_dist <- distance_stan_data$max_dist / 100

# get rid of centre/scale attributes for modelling
distance_stan_data$pitch <- distance_stan_data$pitch[,1]
distance_stan_data$mass <- distance_stan_data$mass[,1]






####### Run Model #################################

inits <- generate_distance_inits(n_chains = n_chains,
                                 napops_skip = c("BITH", "HASP", "KIWA", "LCTH", "LEPC", "SPOW"),
                                 phylo_corr = phylo_corr,
                                 param = "mixed",
                                 species_cp = distance_stan_data$species_cp,
                                 species_ncp = distance_stan_data$species_ncp)

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
