####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 2-generate-distance-cv-folds.R
# Created October 2023
# Last Updated October 2023

####### Import Libraries and External Files #######

library(dplyr)

source("src/functions/subset-distance-data.R")

####### Load Data #################################

load("data/generated/distance_stan_data.rda")

####### Set Constants #############################

k <- 5

####### Allocate CV Numbers #######################

#' First we must drop all species which are just not going to be involved
#' in the cross validation at all. This will be the species for which we are
#' generated predictions
pred_drops <- c("LCTH", "LEPC", "HASP", "SPOW", "KIWA", "BITH")
dis_data <- subset_distance_data(distance_stan_data = distance_stan_data,
                                 sps = setdiff(unique(distance_stan_data$sp_list),
                                               pred_drops))
rm(distance_stan_data) #clear up some mem

#' Create a dataframe that is just the species (for each observation) and
#' fold number
cv_folds <- data.frame(Species = dis_data$sp_list,
                      cv_fold = NA)

for (sp in unique(cv_folds$Species))
{
  n_obs <- nrow(cv_folds[which(cv_folds$Species == sp), ])
  
  obs_per_fold <- ceiling(n_obs / k)
  
  fold <- NULL
  for (i in 1:k)
  {
    fold <- c(fold, rep(i, obs_per_fold))
  }
  fold <- sample(fold[1:n_obs])
  
  cv_folds[which(cv_folds$Species == sp), "cv_fold"] <- fold
}

####### Output ####################################

write.csv(cv_folds, file = "data/generated/distance_cv_folds.csv", row.names = FALSE)
