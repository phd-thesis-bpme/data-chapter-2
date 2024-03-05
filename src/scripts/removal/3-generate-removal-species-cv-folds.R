####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 3-generate-removal-species-cv-folds.R
# Created March 2024
# Last Updated March 2024

####### Import Libraries and External Files #######

library(ape)
library(phytools)

####### Load Data #################################

load("data/generated/removal_stan_data_cv.rda")
phylo_tree <- ape::read.nexus(file = "data/raw/all_species.nex")
traits <- read.csv("data/raw/traits.csv")

####### Set Constants #############################

k <- 5

####### Allocate CV Numbers #######################

tree <- consensus.edges(phylo_tree, method = "mean.edge") |>
keep.tip(dimnames(removal_stan_data_cv$phylo_corr)[[1]])

cut <- cutree((tree), k = k)

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
