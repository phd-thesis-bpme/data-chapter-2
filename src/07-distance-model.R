####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 07-distance-model.R
# Created October 2022
# Last Updated May 2023

####### Import Libraries and External Files #######

library(cmdstanr)
library(dplyr)

####### Load Data #################################

load("data/generated/distance_stan_data.rda")

####### Set Constants #############################

distance_stan_data$grainsize <- 1
distance_stan_data$lambda <- 0

# Scale the maximum distances to units of KM for computational ease
#distance_stan_data$max_dist <- distance_stan_data$max_dist / 1000

# Stan settings
n_iter <- 100
n_warmup <- 250
n_chains <- 1
refresh <- 10
threads_per_chain <- 16

####### Subset Data for Testing ###################

orig_data_df <- data.frame(species = distance_stan_data$sp_list[,1],
                           sp_n = distance_stan_data$species)
### NOTE there are only 322 species in this list.
### is n_species wrong (== 323)
sp_df <- orig_data_df %>%
  group_by(species,sp_n) %>%
  summarise(n_counts = n())

sub_data <- data.frame(mig_strat = distance_stan_data$mig_strat,
                       habitat = distance_stan_data$habitat,
                       pitch = distance_stan_data$pitch,
                       mass = distance_stan_data$mass,
                       sp_n = 1:length(distance_stan_data$mig_strat))
sp_df <- left_join(sp_df,sub_data)

 sps <- c(
#          "WOTH",
#         # "BCCH",
          "WCSP",
        "FLSJ",
         "KIWA",
         "HASP",
"SPOW",
"MEJA",
"BOSP",
        # "BAOR",
        # "EAME",
        # "HAWO",
         "GCFL"
         #"BOBO"
        )
sp_sel <- sp_df[which(sp_df$species %in% sps),]
sp_sel[,"new_sp_n"] <- 1:length(sps) # new species indicators

#identify which of the original data objects should be retained
data_sel <- which(orig_data_df$sp_n %in% sp_sel$sp_n)
sp_data_repl <- orig_data_df %>%
  left_join(.,sp_sel,
            by = "sp_n")

distance_stan_data2 <- distance_stan_data; rm(distance_stan_data); gc();
distance_stan_data2$bands_per_sample <- distance_stan_data2$bands_per_sample[data_sel]
distance_stan_data2$abund_per_band <- distance_stan_data2$abund_per_band[data_sel,]
distance_stan_data2$max_dist <- distance_stan_data2$max_dist[data_sel,]
distance_stan_data2$species <- sp_data_repl[data_sel,"new_sp_n"]
distance_stan_data2$mig_strat <- as.integer(unlist(sp_sel[,"mig_strat"]))
distance_stan_data2$habitat <- as.integer(unlist(sp_sel[,"habitat"]))
distance_stan_data2$mass <- as.numeric(unlist(sp_sel[,"mass"]))
distance_stan_data2$pitch <- as.numeric(unlist(sp_sel[,"pitch"]))

distance_stan_data2$grainsize <- 1
distance_stan_data2$n_species <- length(sps)
distance_stan_data2$n_samples <- length(data_sel)

distance_stan_data2$sp_list <- NULL
distance_stan_data2$phylo_corr <- NULL
distance_stan_data2$max_dist <- distance_stan_data2$max_dist / 100

# # test area for mixture of centred and non-centred parameters
distance_stan_data2$n_species_ncp <- 6
distance_stan_data2$n_species_cp <- 2
distance_stan_data2$species_ncp <- c(1,2,4:7)
distance_stan_data2$species_cp <- c(3,8)

####### Run Model #################################

# # get rid of centre/scale attributes for modelling
# distance_stan_data2$pitch <- NULL#distance_stan_data$pitch[,1]
# distance_stan_data2$mass <- NULL#distance_stan_data$mass[,1]
# distance_stan_data2$n_mig_strat <- NULL
# distance_stan_data2$mig_strat <- NULL
# distance_stan_data2$n_habitat <- NULL
# distance_stan_data2$habitat <- NULL

model_file <- cmdstan_model(stan_file = "models/distance_mixed.stan",
                            cpp_options = list(stan_threads = TRUE))

stan_run <- model_file$sample(
  data = distance_stan_data2,
  iter_warmup = n_warmup,
  iter_sampling = n_iter,
  chains = n_chains,
  parallel_chains = n_chains,
  refresh = refresh,
  threads_per_chain = threads_per_chain,
  adapt_delta = 0.95
)
stan_run$save_object(file = paste0("output/model_runs/distance_predictions.RDS"))

