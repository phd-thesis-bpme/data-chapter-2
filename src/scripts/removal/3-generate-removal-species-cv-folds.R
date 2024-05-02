####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 3-generate-removal-species-cv-folds.R
# Created March 2024
# Last Updated March 2024

####### Import Libraries and External Files #######

library(ape)
library(dplyr)

####### Load Data #################################

load("data/generated/removal_stan_data_cv.rda")
phylo_tree <- ape::read.nexus(file = "data/raw/all_species.nex")
traits <- read.csv("data/raw/traits.csv")
binomial_names <- read.csv("data/generated/binomial_names.csv")

####### Set Constants #############################

k <- 10

# Shamelessly stolen from https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

####### Allocate CV Numbers #######################

sp_to_run <- rownames(removal_stan_data_cv$phylo_corr)

species_groups_df <- data.frame(Species = sp_to_run)

for (i in 1:length(phylo_tree))
{
  t <- cutree(as.hclust.phylo(keep.tip(phylo_tree[[i]], sp_to_run)), k = k)
  sp_names <- names(t)
  temp <- data.frame(Species = sp_names, tree = unname(t))
  names(temp)[2] <- paste0("tree_", i)
  species_groups_df <- merge(species_groups_df, temp, by = "Species")
}

species_groups_df$Consensus <- NA

for (i in 1:nrow(species_groups_df))
{
  species_groups_df$Consensus[i] <- Mode(as.vector(unname(species_groups_df[i, 2:2001])))[[1]]
}

species_groups_df <- species_groups_df[, c("Species", "Consensus")]
species_groups_df$Species <- gsub("_", " ", species_groups_df$Species)
species_groups_df <- merge(species_groups_df, binomial_names[, c("Scientific_BT", "Code")],
                           by.x = "Species", by.y = "Scientific_BT")
species_groups_df <- merge(species_groups_df, traits[, c("Code", "Migrant")],
                           by = "Code")
species_groups_df$Group <- paste0(species_groups_df$Consensus, "-", species_groups_df$Migrant)

species_groups_df$cv_fold <- NA
for (g in unique(species_groups_df$Group))
{
  n_obs <- nrow(species_groups_df[which(species_groups_df$Group == g), ])
  
  obs_per_fold <- ceiling(n_obs / k)
  
  fold <- NULL
  for (f in 1:k)
  {
    fold <- c(fold, rep(f, obs_per_fold))
  }
  if (n_obs <= 20)
  {
    fold <- sample(fold, size = n_obs)
  }else {
    fold <- sample(fold[1:n_obs])
  }
  
  species_groups_df[which(species_groups_df$Group == g), "cv_fold"] <- fold
}

cv_folds <- data.frame(Code = removal_stan_data_cv$sp_list$Species)
cv_folds <- dplyr::left_join(cv_folds, species_groups_df[, c("Code", "cv_fold")], 
        by = "Code")
names(cv_folds) <- c("Species", "cv_fold")

####### Output ####################################

write.csv(cv_folds, file = "data/generated/removal_cv_folds_species.csv", row.names = FALSE)
