####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 03-distance-model.R
# Created April 2022
# Last Updated April 2022

####### Import Libraries and External Files #######

library(plyr)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

####### Read Data #################################

load("data/dist_count_matrix.rda")
load("data/dist_design.rda")
load("output/turdidae_corr_matrix.rda")
binomial <- read.csv("output/binomial_names.csv")

####### Wrangle Data for Modelling ################

# Most of this code adopted from Edwards et al. 2022

max_bands <- ncol(dist_design) - 2
count_names <- c("Sample_ID", "Species", "Distance_Method",
                 paste0(rep("Int", times = max_bands), 1:max_bands))
names(dist_count_matrix) <- count_names

design_names <- c("Distance_Method", "Max_Distance",
                  paste0(rep("Interval", times = max_bands), 1:max_bands))
names(dist_design) <- design_names

# Join data
count_design <- plyr::join(dist_count_matrix, dist_design,
                           by = "Distance_Method", type = "left")

# Filter out methods with only 1 interval
to_remove <- which(is.na(count_design$Interval2))
count_design <- count_design[-c(to_remove), ]

# Create separate data frames
counts <- count_design[,c("Sample_ID", "Species", "Distance_Method", count_names[4:length(count_names)])]
design <- count_design[,c("Sample_ID", "Species", "Distance_Method", design_names[3:length(design_names)])]
names(design) <- names(counts)
col_names <- names(counts)[4:length(names(counts))]

# Change 0s to NA in counts table where appropriate based on design table
for (i in col_names)
{
  indices <- which(counts[, i] == 0)
  
  counts[indices, i] <- ifelse(is.na(design[indices, i]), NA, 0)
}

# Re-order the species so that they match the correlation matrix
species <- sort(as.character(unique(counts$Species)))
bin_turdidae <- binomial[which(binomial$Code %in% species), ]
bin_turdidae$Scientific <- gsub(" ", "_", bin_turdidae$Scientific)
bin_turdidae <- bin_turdidae[match(colnames(corr_matrix), bin_turdidae$Scientific),]
species <- bin_turdidae$Code

# List of input counts
Y <- vector(mode = "list", length = length(species))
names(Y) <- species

# List of input design
D <- vector(mode = "list", length = length(species))
names(D) <- species

# Species indicator list
sp_list <- vector(mode = "list", length = length(species))
names(sp_list) <- species

for (s in species)
{
  n_counts <- nrow(counts[counts$Species == s, ])
  
  Y[[s]] <- as.matrix(counts[counts$Species==s, col_names])
  D[[s]] <- as.matrix(design[design$Species==s, col_names])
  sp_list[[s]] <- data.frame(Species = rep(s, n_counts))
}

Y <- Y[lengths(Y) != 0]; Y <- do.call(rbind, Y)
D <- D[lengths(D) != 0]; D <- do.call(rbind, D)
sp_list <- sp_list[lengths(sp_list) != 0]; sp_list <- do.call(rbind, sp_list)

#' Corresponds with "bands_per_sample" in distance.stan
dist_bands_per_sample <- unname(apply(Y, 1, function(x) sum(!is.na(x))))

#' Total species abundance per sampling event.
#' I.e., this is the sum of Y_sik over k
#' Corresponds with "abund_per_sample" in distance.stan
total_abund_per_sample <- unname(apply(Y, 1, function(x) sum(x, na.rm = TRUE)))

#' Factored list of species
#' Corresponds with "species" in distance.stan
sp_list_numeric <- as.numeric(factor(sp_list[,1],
                                     levels = species))

#' Corresponds with "abund_per_band" in distance.stan
abundance_per_band <- Y
abundance_per_band[is.na(abundance_per_band)] <- 0

#' Corresponds with "max_dist" in distance.stan
max_dist <- D
max_dist[is.na(max_dist)] <- 0

n_samples <- nrow(Y)
n_species <- length(unique(sp_list[,1]))
max_intervals <- ncol(Y)

stan_data <- list(n_samples = n_samples,
                  n_species = n_species,
                  max_intervals = max_intervals,
                  species = sp_list_numeric,
                  abund_per_band = abundance_per_band,
                  abund_per_sample = total_abund_per_sample,
                  bands_per_sample = dist_bands_per_sample,
                  max_dist = max_dist,
                  phylo_corr = corr_matrix)

####### Output ####################################

model <- stan_model(file = "models/distance.stan")

stan_job <- sampling(model,
                     data = stan_data,
                     verbose = TRUE,
                     chains = 4,
                     iter = 1000,
                     warmup = 500,
                     cores = 4,
                     pars = c("log_tau", "sigma", "mu"),
                     control = list(adapt_delta = 0.8,
                                    max_treedepth = 15))

####### Output ####################################

save(stan_job, file = "output/distance_turdidae_model.rda")
save(stan_data, file = "output/distance_turdidae_data.rda")