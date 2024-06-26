####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# posthoc/3-kfold-cv-figure.R
# Created December 2023
# Last Updated June 2024

####### Import Libraries and External Files #######

library(cmdstanr)
library(bayesplot)
library(ggplot2)
library(ggpubr)
library(ape)
library(ggtree)
theme_set(theme_pubclean())
bayesplot::color_scheme_set("red")

source("src/functions/rem-lik.R")
source("src/functions/dis-lik.R")

####### Read Data #################################

load("data/generated/removal_stan_data_cv.rda")
load("data/generated/distance_stan_data_cv.rda")
cv_folds_rem <- read.csv("data/generated/removal_cv_folds.csv")
cv_folds_dis <- read.csv("data/generated/distance_cv_folds.csv")
phylo_tree <- ape::read.nexus(file = "data/raw/all_species.nex")
binomial <- read.csv("data/generated/binomial_names.csv")
traits <- read.csv("data/raw/traits.csv")

####### Removal Model #############################

#' Create dataframe that will connect species, data point,
#' and lppd
lppd <- data.frame(Species = cv_folds_rem$Species,
                   CV_Fold = cv_folds_rem$cv_fold,
                   SS_LPPD = NA,
                   MS_LPPD = NA,
                   Difference = NA)

for (f in 1:max(cv_folds_rem$cv_fold)){
  
  # Extract dataframe of phi from single species and multi species model
  ss_mod <- readRDS(paste0("output/model_runs/cv_removal/ms-vs-ss/ss_fold_",
                           f,
                           ".RDS"))
  ms_mod <- readRDS(paste0("output/model_runs/cv_removal/ms-vs-ss/ms_fold_",
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
                           Mean_Difference = NA,
                           StdErr = NA,
                           Model_Preference = NA)

sp_sample_size <- data.frame(table(removal_stan_data_cv$sp_list$Species))
names(sp_sample_size) <- c("Species", "N")
lppd_summary <- dplyr::left_join(lppd_summary, sp_sample_size, by = "Species")

for (i in 1:nrow(lppd_summary))
{
  sp <- lppd_summary$Species[i]
  lppd_summary$Mean_Difference[i] <- mean(lppd[which(lppd$Species == sp),
                                               "Difference"])
  lppd_summary$StdErr[i] <- sqrt(var(lppd[which(lppd$Species == sp),
                                          "Difference"])) / lppd_summary$N[i]
  
}

lppd_summary[which((lppd_summary$Mean_Difference - lppd_summary$StdErr < 0) &
                     (lppd_summary$Mean_Difference + lppd_summary$StdErr < 0)),
             "Model_Preference"] <- "SS"

lppd_summary[which((lppd_summary$Mean_Difference - lppd_summary$StdErr > 0) &
                     (lppd_summary$Mean_Difference + lppd_summary$StdErr > 0)),
             "Model_Preference"] <- "MS"

preference <- data.frame(table(lppd_summary$Model_Preference))

(lppd_difference_plot_rem <- ggplot(data = lppd_summary, aes(x = Species, y = Mean_Difference)) +
    geom_errorbar(aes(ymin = Mean_Difference - StdErr, ymax = Mean_Difference + StdErr)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    #coord_flip() +
    NULL)

(lppd_diff_vs_n_rem <- ggplot(data = lppd_summary, aes(x = log(N), y = Mean_Difference)) +
    geom_point(aes(color = Model_Preference)) +
    geom_hline(yintercept = 0) +
    xlab("Log(Sample Size)") +
    ylab("Mean LPPD Difference") +
    theme(legend.position = "none") +
    NULL)

lppd_summary <- merge(lppd_summary, binomial[,c("Code", "Scientific_BT")],
                      by.x = "Species", by.y = "Code")
lppd_summary$Scientific_BT <- gsub(" ", "_", lppd_summary$Scientific_BT)
lppd_summary <- merge(lppd_summary, traits[, c("Code", "Migrant")],
                      by.x = "Species", by.y = "Code")
tree <- ape::consensus(phylo_tree)
drops <- binomial[-which(binomial$Code %in% removal_stan_data_cv$sp_list$Species), "Scientific_BT"]
drops <- gsub(" ", "_", drops)
tree <- drop.tip(tree, drops)

(model_pref_tree_rem <- ggtree(tree, branch.length="none") %<+% lppd_summary[, c("Scientific_BT", "Model_Preference", "N", "Migrant")] +
    geom_tiplab(offset = .6, hjust = .5) +
    geom_tippoint(aes(color = Model_Preference, shape = as.factor(Migrant), size = log(N))))

lppd_removal <- lppd
lppd_summary_removal <- lppd_summary

####### Distance Model ############################

#' Create dataframe that will connect species, data point,
#' and lppd
lppd <- data.frame(Species = cv_folds_dis$Species,
                   CV_Fold = cv_folds_dis$cv_fold,
                   SS_LPPD = NA,
                   MS_LPPD = NA,
                   Difference = NA)

for (f in 1:max(cv_folds_dis$cv_fold))
{
  # Extract dataframe of tau from single species and multi species model
  ss_mod <- readRDS(paste0("output/model_runs/cv_distance/ms-vs-ss/ss_fold_",
                           f,
                           ".RDS"))
  ms_mod <- readRDS(paste0("output/model_runs/cv_distance/ms-vs-ss/ms_fold_",
                           f,
                           ".RDS"))
  
  ss_tau_summary <- ss_mod$summary("log_tau")
  ms_tau_summary <- ms_mod$summary("log_tau")
  
  # Subset the data to just the test set
  data_indices <- which(lppd$CV_Fold == f)
  abund <- distance_stan_data_cv$abund_per_band[data_indices, ]
  bands <- distance_stan_data_cv$bands_per_sample[data_indices]
  max_dist <- distance_stan_data_cv$max_dist[data_indices, ] / 100
  species <- distance_stan_data_cv$species[data_indices]
  max_intervals <- distance_stan_data_cv$max_intervals
  
  tau_ss <- exp(ss_tau_summary$mean[species])
  tau_ms <- exp(ms_tau_summary$mean[species])
  
  # Check to see if we can drop any training columns (i.e. if bands per sample is smaller)
  if (max(bands) < max_intervals)
  {
    max_intervals <- max(bands)
    abund <- abund[, 1:max_intervals]
    max_dist <- max_dist[, 1:max_intervals]
  }
  
  lp_ss <- dis_lik(abund = abund,
                   bands = bands,
                   max_dist = max_dist,
                   max_intervals = max_intervals,
                   tau = tau_ss)
  
  lp_ms <- dis_lik(abund = abund,
                   bands = bands,
                   max_dist = max_dist,
                   max_intervals = max_intervals,
                   tau = tau_ms)
  
  lppd$SS_LPPD[data_indices] <- lp_ss
  lppd$MS_LPPD[data_indices] <- lp_ms
}

lppd$Difference <- lppd$MS_LPPD - lppd$SS_LPPD

lppd_summary <- data.frame(Species = unique(lppd$Species),
                           Mean_Difference = NA,
                           StdErr = NA,
                           Model_Preference = NA)

sp_sample_size <- data.frame(table(distance_stan_data_cv$sp_list))
names(sp_sample_size) <- c("Species", "N")
lppd_summary <- dplyr::left_join(lppd_summary, sp_sample_size, by = "Species")

for (i in 1:nrow(lppd_summary))
{
  sp <- lppd_summary$Species[i]
  lppd_summary$Mean_Difference[i] <- mean(lppd[which(lppd$Species == sp),
                                               "Difference"])
  lppd_summary$StdErr[i] <- sqrt(var(lppd[which(lppd$Species == sp),
                                          "Difference"])) / lppd_summary$N[i]
  
}

lppd_summary[which((lppd_summary$Mean_Difference - lppd_summary$StdErr < 0) &
                     (lppd_summary$Mean_Difference + lppd_summary$StdErr < 0)),
             "Model_Preference"] <- "SS"

lppd_summary[which((lppd_summary$Mean_Difference - lppd_summary$StdErr > 0) &
                     (lppd_summary$Mean_Difference + lppd_summary$StdErr > 0)),
             "Model_Preference"] <- "MS"

preference <- data.frame(table(lppd_summary$Model_Preference))

(lppd_difference_plot_dis <- ggplot(data = lppd_summary, aes(x = Species, y = Mean_Difference)) +
    geom_errorbar(aes(ymin = Mean_Difference - StdErr, ymax = Mean_Difference + StdErr)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    #coord_flip() +
    NULL)

(lppd_diff_vs_n_dis <- ggplot(data = lppd_summary, aes(x = log(N), y = Mean_Difference)) +
    geom_point(aes(color = Model_Preference)) +
    geom_hline(yintercept = 0) +
    xlab("Log(Sample Size)") +
    ylab("Mean LPPD Difference") +
    theme(legend.position = "none") +
    NULL
)

####### Output ####################################

write.table(lppd_removal, file = "data/generated/lppd_removal.csv", sep = ",", row.names = FALSE)
write.table(lppd_summary_removal, file = "data/generated/lppd_summary_removal.csv", sep = ",", row.names = FALSE)
write.table(lppd, file = "data/generated/lppd_distance.csv", sep = ",", row.names = FALSE)
write.table(lppd_summary, file = "data/generated/lppd_summary_distance.csv", sep = ",", row.names = FALSE)

tiff(filename = "output/plots/kfold_cv_plot.tiff",
    width = 6, height = 6, units = "in", res = 600)
ggarrange(lppd_diff_vs_n_rem, lppd_diff_vs_n_dis, nrow = 2, labels = c("A", "B"))
dev.off()

png(filename = "output/plots/kfold_cv_plot.png",
     width = 6, height = 6, units = "in", res = 600)
ggarrange(lppd_diff_vs_n_rem, lppd_diff_vs_n_dis, nrow = 2, labels = c("A", "B"))
dev.off()

