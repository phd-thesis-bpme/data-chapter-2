####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 3-generate-distance-species-cv-folds.R
# Created March 2024
# Last Updated April 2024

####### Import Libraries and External Files #######

library(dplyr)

source("src/functions/subset-distance-data.R")

####### Load Data #################################

load("data/generated/distance_stan_data.rda")

####### Set Constants #############################

k <- 10

####### Allocate CV Numbers #######################

#' First we must drop all species which are just not going to be involved
#' in the cross validation at all. This will be the species for which we are
#' generated predictions
pred_drops <- c("LCTH", "LEPC", "HASP", "SPOW", "KIWA", "BITH")
dis_data <- subset_distance_data(distance_stan_data = distance_stan_data,
                                 sps = setdiff(unique(distance_stan_data$sp_list),
                                               pred_drops))
rm(distance_stan_data) #clear up some mem

species_groups_df <- data.frame(Species = setdiff(dis_data$sp_all, pred_drops),
                                Migrant = dis_data$mig_strat,
                                Habitat = dis_data$habitat)
species_groups_df$group <- paste0(species_groups_df$Migrant, "-", species_groups_df$Habitat)
species_groups_df$cv_fold <- NA

for (g in unique(species_groups_df$group))
{
  n_obs <- nrow(species_groups_df[which(species_groups_df$group == g), ])
  
  obs_per_fold <- ceiling(n_obs / k)
  
  fold <- NULL
  for (f in 1:k)
  {
    fold <- c(fold, rep(f, obs_per_fold))
  }
  
  fold <- sample(fold, n_obs)
  
  species_groups_df[which(species_groups_df$group == g), "cv_fold"] <- fold
}

cv_folds <- data.frame(Species = dis_data$sp_list)
cv_folds <- dplyr::left_join(cv_folds, species_groups_df[, c("Species", "cv_fold")],
                             by = "Species")

####### Output ####################################

write.csv(cv_folds, file = "data/generated/distance_cv_folds_species.csv", row.names = FALSE)
