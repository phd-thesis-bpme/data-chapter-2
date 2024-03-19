####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 3-generate-removal-species-cv-folds.R
# Created March 2024
# Last Updated March 2024

####### Import Libraries and External Files #######

library(ape)

####### Load Data #################################

load("data/generated/removal_stan_data_cv.rda")
phylo_tree <- ape::read.nexus(file = "data/raw/all_species.nex")
traits <- read.csv("data/raw/traits.csv")
binomial_names <- read.csv("data/generated/binomial_names.csv")

####### Set Constants #############################

k <- 5

# Shamelessly stolen from https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

####### Allocate CV Numbers #######################





species_groups_df <- data.frame(Species = gsub(" ", "_", x = binomial_names$Scientific_BT))

for (i in 1:length(phylo_tree))
{
  t <- cutree(as.hclust.phylo(phylo_tree[[i]]), k = k)
  sp_names <- names(t)
  temp <- data.frame(Species = sp_names, tree = unname(t))
  names(temp)[2] <- paste0("tree_", i)
  species_groups_df <- merge(species_groups_df, temp, by = "Species")
}

species_groups_df$Consensus <- NA

for (i in 1:nrow(species_groups_df))
{
  species_groups_df$Consensus[i] <- Mode(as.vector(unname(species_groups_df[i, 2:2001])))
}

species_groups_df <- species_groups_df[, c("Species", "Consensus")]











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
