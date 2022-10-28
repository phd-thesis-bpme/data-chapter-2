####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 02-prepare-removal-data.R
# Created August 2022
# Last Updated October 2022

####### Import Libraries and External Files #######

library(ape)
library(Matrix)
library(plyr)
library(magrittr)

####### Read Data #################################

load("data/raw/time_count_matrix.rda")
load("data/raw/time_design.rda")
cv_tree <- load(file = "data/generated/corr_matrix_cv.rda")
pred_tree <- load(file = "data/generated/corr_matrix_predict.rda")
binomial <- read.csv("data/generated/binomial_names.csv")

####### Wrangle Data for Modelling ################

# Most of this code adopted from Edwards et al. 2022

# Drop the 99 column
time_count_matrix <- time_count_matrix[, c(1:ncol(time_count_matrix) - 1)]

max_bands <- ncol(time_design) - 2
count_names <- c("Sample_ID", "Species", "Time_Method",
                 paste0(rep("Int", times = max_bands), 1:max_bands))
names(time_count_matrix) <- count_names

design_names <- c("Time_Method", "Max_Duration",
                  paste0(rep("Interval", times = max_bands), 1:max_bands))
names(time_design) <- design_names

# Join data
count_design <- plyr::join(time_count_matrix, time_design,
                           by = "Time_Method", type = "left")

# Filter out methods with only 1 interval
to_remove <- which(is.na(count_design$Interval2))
count_design <- count_design[-c(to_remove), ]

# Create separate data frames
counts <- count_design[,c("Sample_ID", "Species", "Time_Method", count_names[4:length(count_names)])]
design <- count_design[,c("Sample_ID", "Species", "Time_Method", design_names[3:length(design_names)])]
names(design) <- names(counts)
col_names <- names(counts)[4:length(names(counts))]

# Change 0s to NA in counts table where appropriate based on design table
for (i in col_names)
{
  indices <- which(counts[, i] == 0)
  
  counts[indices, i] <- ifelse(is.na(design[indices, i]), NA, 0)
}

# Get species codes and order of species codes for the Cross-validation matrix
species_cv <- gsub(pattern = "_", 
                   replacement = " ", 
                   x = rownames(corr_matrix_cv))
binomial_cv <- binomial[which(binomial$Scientific_BT %in% species_cv), ]
binomial_cv <- binomial_cv[match(species_cv, binomial_cv$Scientific_BT), ]
species_cv_code <- binomial_cv$Code

# Get species codes and order of species codes for the Prediction Matrix
species_pred <- gsub(pattern = "_",
                     replacement= " ",
                     x = rownames(corr_matrix_predict))
binomial_pred <- binomial[which(binomial$Scientific_BT %in% species_pred), ]
binomial_pred <- binomial_pred[match(species_pred, binomial_pred$Scientific_BT), ]
species_pred_code <- binomial_pred$Code

#' All of the following values will have values specific to cross validation list
#' and the prediction list. So those will be noted e.g. Y_cv or Y_pred

# List of input counts
Y_cv <- vector(mode = "list", length = length(species_cv_code))
names(Y_cv) <- species_cv_code

Y_pred <- vector(mode = "list", length = length(species_pred_code))
names(Y_pred) <- species_pred_code

# List of input design
D_cv <- vector(mode = "list", length = length(species_cv_code))
names(D_cv) <- species_cv_code

D_pred <- vector(mode = "list", length = length(species_pred_code))
names(D_pred) <- species_pred_code

# Species indicator list
sp_list_cv <- vector(mode = "list", length = length(species_cv_code))
names(sp_list_cv) <- species_cv_code

sp_list_pred <- vector(mode = "list", length = length(species_pred_code))
names(sp_list_pred) <- species_pred_code

for (s in species_cv_code)
{
  n_counts <- nrow(counts[counts$Species == s, ])
  
  Y_cv[[s]] <- as.matrix(counts[counts$Species==s, col_names])
  D_cv[[s]] <- as.matrix(design[design$Species==s, col_names])
  sp_list_cv[[s]] <- data.frame(Species = rep(s, n_counts))
}

for (s in species_pred_code)
{
  n_counts <- nrow(counts[counts$Species == s, ])
  
  Y_pred[[s]] <- as.matrix(counts[counts$Species==s, col_names])
  D_pred[[s]] <- as.matrix(design[design$Species==s, col_names])
  sp_list_pred[[s]] <- data.frame(Species = rep(s, n_counts))
}

Y_cv <- Y_cv[lengths(Y_cv) != 0]; Y_cv <- do.call(rbind, Y_cv)
Y_pred <- Y_pred[lengths(Y_pred) != 0]; Y_pred <- do.call(rbind, Y_pred)

D_cv <- D_cv[lengths(D_cv) != 0]; D_cv <- do.call(rbind, D_cv)
D_pred <- D_pred[lengths(D_pred) != 0]; D_pred <- do.call(rbind, D_pred)

sp_list_cv <- do.call(rbind, sp_list_cv)
sp_list_pred <- do.call(rbind, sp_list_pred)

#' Corresponds with "bands_per_sample" in removal.stan
time_bands_per_sample_cv <- unname(apply(Y_cv, 1, function(x) sum(!is.na(x))))
time_bands_per_sample_pred <- unname(apply(Y_pred, 1, function(x) sum(!is.na(x))))

#' Total species abundance per sampling event.
#' I.e., this is the sum of Y_sij over j
#' Corresponds with "abund_per_sample" in removal.stan
total_abund_per_sample_cv <- unname(apply(Y_cv, 1, function(x) sum(x, na.rm = TRUE)))
total_abund_per_sample_pred <- unname(apply(Y_pred, 1, function(x) sum(x, na.rm = TRUE)))

#' Factored list of species
#' Corresponds with "species" in removal.stan
sp_cv_numeric <- data.frame(Species = species_cv_code,
                            num = seq(1, length(species_cv_code)))
sp_cv_numeric <- join(sp_list_cv, sp_cv_numeric, by= "Species")

sp_pred_numeric <- data.frame(Species = species_pred_code,
                            num = seq(1, length(species_pred_code)))
sp_pred_numeric <- join(sp_list_pred, sp_pred_numeric, by= "Species")

#' Corresponds with "abund_per_band" in removal.stan
abundance_per_band_cv <- Y_cv
abundance_per_band_cv[is.na(abundance_per_band_cv)] <- 0

abundance_per_band_pred <- Y_pred
abundance_per_band_pred[is.na(abundance_per_band_pred)] <- 0

#' Corresponds with "max_time" in removal.stan
max_time_cv <- D_cv
max_time_cv[is.na(max_time_cv)] <- 0

max_time_pred <- D_pred
max_time_pred[is.na(max_time_pred)] <- 0

n_samples_cv <- nrow(Y_cv)
n_samples_pred <- nrow(Y_pred)

n_species_cv <- nrow(binomial_cv)
n_species_pred <- nrow(binomial_pred)

max_intervals_cv <- ncol(Y_cv)
max_intervals_pred <- ncol(Y_pred)

removal_stan_data_cv <- list(n_samples = n_samples_cv,
                          n_species = n_species_cv,
                          max_intervals = max_intervals_cv,
                          species = sp_cv_numeric$num,
                          abund_per_band = abundance_per_band_cv,
                          bands_per_sample = time_bands_per_sample_cv,
                          max_time = max_time_cv,
                          phylo_corr = cv_tree)

removal_stan_data_pred <- list(n_samples = n_samples_pred,
                             n_species = n_species_pred,
                             max_intervals = max_intervals_pred,
                             species = sp_pred_numeric$num,
                             abund_per_band = abundance_per_band_pred,
                             bands_per_sample = time_bands_per_sample_pred,
                             max_time = max_time_pred,
                             phylo_corr = pred_tree)

####### Output ####################################
save(removal_stan_data_cv, file = "data/generated/removal_stan_data_cv.rda")
save(removal_stan_data_pred, file = "data/generated/removal_stan_data_pred.rda")

