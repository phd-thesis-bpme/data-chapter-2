####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 3-generate-distance-species-cv-folds.R
# Created March 2024
# Last Updated April 2024

####### Load Data #################################

load("data/generated/distance_stan_data_cv.rda")

####### Set Constants #############################

k <- 10

####### Allocate CV Numbers #######################

species_groups_df <- data.frame(Species = distance_stan_data_cv$sp_all, 
                                Migrant = distance_stan_data_cv$mig_strat,
                                Habitat = distance_stan_data_cv$habitat)
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

cv_folds <- data.frame(Species = distance_stan_data_cv$sp_list)
cv_folds <- dplyr::left_join(cv_folds, species_groups_df[, c("Species", "cv_fold")],
                             by = "Species")

####### Output ####################################

write.csv(cv_folds, file = "data/generated/distance_cv_folds_species.csv", row.names = FALSE)
