####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 07-generate-distance-cv-folds.R
# Created October 2023
# Last Updated October 2023

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

# Create a dataframe of traits from the current distance data, including species name
traits <- data.frame(Mig = dis_data$mig_strat,
                     Habitat = dis_data$habitat,
                     Mass = dis_data$mass,
                     Pitch = dis_data$pitch,
                     Species = unique(dis_data$sp_list))

# Create an indicator for each Mig x Hab association
traits$Mig_Hab <- paste0(traits$Mig, "-", traits$Habitat)
counts_per_group <- data.frame(table(traits$Mig_Hab))

traits$cv_fold <- NA

for (i in 1:nrow(counts_per_group))
{
  group <- counts_per_group$Var1[i]
  n_sp <- counts_per_group$Freq[i]
  
  cv_allocs <- rep(seq(1,k), times = ceiling(n_sp / k))
  cv_allocs <- cv_allocs[1:n_sp]
  cv_alloc_random <- sample(cv_allocs)
  
  traits[which(traits$Mig_Hab == group), "cv_fold"] <- cv_alloc_random
}

cv_folds <- traits[, c("Species", "cv_fold")]

####### Output ####################################

write.csv(cv_folds, file = "data/generated/distance_cv_folds.csv", row.names = FALSE)
