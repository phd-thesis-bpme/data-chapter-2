####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# removal-cv-analysis.R
# Created October 2023
# Last Updated October 2023

####### Import Libraries and External Files #######

library(cmdstanr)
library(bayesplot)
library(ggplot2)
library(ggpubr)
theme_set(theme_pubclean())
bayesplot::color_scheme_set("red")

source("src/functions/rem-lik.R")

####### Read Data #################################

load("data/generated/removal_stan_data_cv.rda")
cv_folds <- read.csv("data/generated/removal_cv_folds.csv")

####### Main Code #################################

#' Create dataframe that will connect species, data point,
#' and lppd
lppd <- data.frame(Species = cv_folds$Species,
                   CV_Fold = cv_folds$cv_fold,
                   SS_LPPD = NA,
                   MS_LPPD = NA,
                   Difference = NA)

for (f in 1:max(cv_folds$cv_fold))
{
  # Extract dataframe of phi from single species and multi species model
  ss_mod <- readRDS(paste0("output/model_runs/cv_removal/ss_fold_",
                    f,
                    ".RDS"))
  ms_mod <- readRDS(paste0("output/model_runs/cv_removal/ms_fold_",
                           f,
                           ".RDS"))
  
  ss_phi_summary <- ss_mod$summary("log_phi")
  ms_phi_summary <- ms_mod$summary("log_phi")
  
  # Subset the data to just the test set
  data_indices <- which(lppd$CV_Fold == f)
  abund <- removal_stan_data_cv$abund_per_band[data_indices, ]
  bands <- removal_stan_data_cv$bands_per_sample[data_indices]
  max_time <- removal_stan_data_cv$max_time[data_indices, ]
  species <- removal_stan_data_cv$species[data_indices]
  max_intervals <- removal_stan_data_cv$max_intervals

  phi_ss <- exp(ss_phi_summary$mean[species])
  phi_ms <- exp(ms_phi_summary$mean[species])
  
  # Check to see if we can drop any training columns (i.e. if bands per sample is smaller)
  if (max(bands) < max_intervals)
  {
    max_intervals <- max(bands)
    abund <- abund[, 1:max_intervals]
    max_time <- max_time[, 1:max_intervals]
  }
  
  lp_ss <- rem_lik(abund = abund,
                   bands = bands,
                   max_time = max_time,
                   max_intervals = max_intervals,
                   phi = phi_ss)
  
  lp_ms <- rem_lik(abund = abund,
                   bands = bands,
                   max_time = max_time,
                   max_intervals = max_intervals,
                   phi = phi_ms)
  
  lppd$SS_LPPD[data_indices] <- lp_ss
  lppd$MS_LPPD[data_indices] <- lp_ms
}

lppd$Difference <- lppd$MS_LPPD - lppd$SS_LPPD

lppd_summary <- data.frame(Species = unique(lppd$Species),
                           Median_LPPD = NA)

sp_sample_size <- data.frame(table(removal_stan_data_cv$sp_list$Species))
names(sp_sample_size) <- c("Species", "N")

lppd_summary <- dplyr::left_join(lppd_summary, sp_sample_size, by = "Species")
for (i in 1:nrow(lppd_summary))
{
  sp <- lppd_summary$Species[i]
  lppd_summary$Median_LPPD[i] <- mean(lppd[which(lppd$Species == sp),
                                             "Difference"])
}

plot(lppd_summary$Median_LPPD ~ lppd_summary$N)




cv_model <- cmdstan_model(stan_file = "models/cv_difference.stan")
cv_model_data <- list(N = nrow(lppd),
                      n_species = length(unique(lppd$Species)),
                      difference = lppd$Difference,
                      species = removal_stan_data_cv$species)

cv_model_run <- cv_model$sample(
  data = cv_model_data,
  iter_warmup = 250,
  iter_sampling = 500,
  chains = 4,
  parallel_chains = 4,
  refresh = 10
)

saveRDS(cv_model_run, file = "output/model_runs/cv_removal/cv_model_removal.RDS")

cv_model_summary <- cv_model_run$summary(c("nu", "overall_difference", "species_difference", "sigma"))

####### Output ####################################

write.table(lppd, file = "data/generated/lppd_removal.csv", sep = ",", row.names = FALSE)
