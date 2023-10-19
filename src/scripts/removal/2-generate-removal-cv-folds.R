####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 2-generate-removal-cv-folds.R
# Created October 2023
# Last Updated October 2023

####### Load Data #################################

load("data/generated/removal_stan_data_cv.rda")

####### Set Constants #############################

k <- 5

####### Allocate CV Numbers #######################

#' Create a dataframe that is just the species (for each observation) and
#' fold number
cv_folds <- data.frame(Species = removal_stan_data_cv$sp_list,
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

write.csv(cv_folds, file = "data/generated/removal_cv_folds.csv", row.names = FALSE)
