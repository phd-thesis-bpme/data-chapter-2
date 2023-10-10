####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 11-distance-cv-analysis.R
# Created October 2023
# Last Updated October 2023
####### Import Libraries and External Files #######

library(cmdstanr)
library(bayesplot)
library(plyr)
library(ggplot2)
library(ggpubr)
theme_set(theme_pubclean())
bayesplot::color_scheme_set("red")

source("src/functions/subset-distance-data.R")
source("src/functions/dis-lik.R")

####### Read Data #################################

cv_folds <- read.csv("data/generated/distance_cv_folds.csv")
load("data/generated/distance_stan_data.rda")

####### Data Wrangling ############################

#' First we must drop all species which are just not going to be involved
#' in the cross validation at all. This will be the species for which we are
#' generated predictions
pred_drops <- c("LCTH", "LEPC", "HASP", "SPOW", "KIWA", "BITH")
dis_data <- subset_distance_data(distance_stan_data = distance_stan_data,
                                 sps = setdiff(unique(distance_stan_data$sp_list),
                                               pred_drops))
rm(distance_stan_data)

cv_folds$Index <- seq(1, nrow(cv_folds))
cv_folds$EDR <- NA

for (i in 1:max(cv_folds$cv_fold))
{
  fold <- readRDS(paste0("output/model_runs/cv_distance/fold_",
                         i,
                         ".RDS"))
  
  # Add the mean estimated EDR to the dataframe above
  log_tau_df <- fold$summary("log_tau")
  fold_indices <- which(cv_folds$cv_fold == i)
  cv_folds[fold_indices, "EDR"] <- exp(log_tau_df[fold_indices, "mean"]) 
  
  # Now calculate the median LPMF for each species
  sp_fold <- cv_folds[which(cv_folds$cv_fold %in% i), c("Species", "Index")]
  
  for (sp_index %in% sp_fold$Index)
  {
    data_indices <- which(dis_data$species == sp_index)
    abund <- dis_data$abund_per_band[data_indices, ]
    bands <- dis_data$bands_per_sample[data_indices]
    max_dist <- dis_data$max_dist[data_indices, ] / 100
    max_intervals <- dis_data$max_intervals
    
    # Check to see if we can drop any training columns (i.e. if bands per sample is smaller)
    if (max(bands) < max_intervals)
    {
      max_intervals <- max(bands)
      abund <- abund[, 1:max_intervals]
      max_dist <- max_dist[, 1:max_intervals]
    }
    
    tau_df <- fold$draws(variables = paste0("log_tau[",
                                            sp_index,
                                            "]"),
                         format = "df")
    
    loglik_vector <- vector(mode = "double", length = nrow(tau_df))
    
    # This can probably be vectorized with apply() or something.
    for (j in 1:length(loglik_vector))
    {
      loglik_vector[j] <- median(dis_lik(abund = abund,
              bands = bands,
              max_dist = max_dist,
              max_intervals = max_intervals,
              tau = exp(tau_df$`log_tau[7]`[j])))
    }
    
    #' TO DO: Need to either find a way to keep track of the loglik_vector for
    #' each species (perhaps as a n_draws X n_species matrix), or just track the
    #' mean/median and 95% intervals of those vectors?
  }
}
