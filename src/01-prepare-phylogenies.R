####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 01-prepare-phylogenies.R
# Created August 2022
# Last Updated October 2022

####### Import Libraries and External Files #######

library(ape)
library(napops)
library(magrittr)

source("src/functions/generate-phylo-corr.R")

####### Read Data #################################

phylo_tree <- ape::read.nexus(file = "data/raw/all_species.nex")
binomial <- read.csv("data/generated/binomial_names.csv")

####### Main Code #################################

#' Phylo distance matrix for cross-validation and precision comparison purposes.
#' I.e., only take species that were a) modelled in Edwards et al. 2022 NA-POPS
#' paper, and b) had sufficient data for BOTH removal and distance
napops_species <- napops::list_species()
napops_species <- napops_species[which(napops_species$Removal == 1 &
                                         napops_species$Distance == 1), ]

# This will drop Woodhouse's Scrub Jay, Pacific Wren, and Sagebrush Sparrow :(
sp_cv <- binomial[which(binomial$Code %in% napops_species$Species), ]
dropped_species <- gsub(" ", "_", setdiff(binomial$Scientific_BT, sp_cv$Scientific_BT))

# Generate the correlation matrix for the cross validation (takes a long time)
corr_matrix_cv <- generate_phylo_corr(phylo_tree = phylo_tree, drops = dropped_species)

#' Now we will create the phylo distance matrix for the sample size relaxation
#' and prediction portion of the experiments. Here, we will add in 5 species to
#' list of species to keep

sp_new <- c("LeConte's Thrasher", "Bicknell's Thrush", 
            "Lesser Prairie-Chicken", "Harris's Sparrow", 
            "Kirtland's Warbler", "Tricolored Blackbird")


sp_predict <- rbind(sp_cv,
                    binomial[which(binomial$English %in% sp_new), ])
dropped_species <- gsub(" ", "_", setdiff(binomial$Scientific_BT, sp_predict$Scientific_BT))

# Generate the correlation matrix for the predictions (takes a long time)
corr_matrix_predict <- generate_phylo_corr(phylo_tree = phylo_tree, drops = dropped_species)

####### Output ####################################

save(corr_matrix_cv, file = "data/generated/corr_matrix_cv.rda")
save(corr_matrix_predict, file = "data/generated/corr_matrix_predict.rda")
