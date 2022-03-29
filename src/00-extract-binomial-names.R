####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 00-extract-binomial-names.R
# Created March 2022
# Last Updated March 2022

####### Import Libraries and External Files #######


####### Read Data #################################

#' This file comes directly from NA-POPS: analysis/01-combine-data.R
#' which is not publicly available. Will be loaded in as "project_counts"
load("data/counts.rda")

codes <- read.csv("util/IBP-Alpha-Codes20.csv")

####### Main Code #################################

# Extract all binomial names for species that NA-POPS has data for
codes_red <- codes[which(codes$SPEC %in% unique(project_counts$Species)), ]
binomial <- codes_red[, c("SPEC", "COMMONNAME", "SCINAME")]

####### Output ####################################

write.table(binomial,
            file = "output/binomial_names.csv",
            sep = ",",
            row.names = FALSE)
