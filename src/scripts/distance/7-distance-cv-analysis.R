####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 7-distance-cv-analysis.R
# Created October 2023
# Last Updated October 2023
####### Import Libraries and External Files #######

library(cmdstanr)
library(bayesplot)
library(plyr)
library(ggplot2)
library(ggpubr)
theme_set(theme_pubclean())
bayesplot::color_scheme_set("red")

source("src/functions/subset-distance-data.R")
source("src/functions/dis-lik.R")

####### Read Data #################################

cv_folds <- read.csv("data/generated/distance_cv_folds.csv")
load("data/generated/distance_stan_data.rda")
traits <- read.csv("data/raw/traits.csv")

####### Data Wrangling ############################

#' First we must drop all species which are just not going to be involved
#' in the cross validation at all. This will be the species for which we are
#' generated predictions
pred_drops <- c("LCTH", "LEPC", "HASP", "SPOW", "KIWA", "BITH")
dis_data <- subset_distance_data(distance_stan_data = distance_stan_data,
                                 sps = setdiff(unique(distance_stan_data$sp_list),
                                               pred_drops))
rm(distance_stan_data)

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
  # Extract dataframe of tau from single species and multi species model
  ss_mod <- readRDS(paste0("output/model_runs/cv_distance/ss_fold_",
                           f,
                           ".RDS"))
  ms_mod <- readRDS(paste0("output/model_runs/cv_distance/ms_fold_",
                           f,
                           ".RDS"))
  
  ss_tau_summary <- ss_mod$summary("log_tau")
  ms_tau_summary <- ms_mod$summary("log_tau")
  
  # Subset the data to just the test set
  data_indices <- which(lppd$CV_Fold == f)
  abund <- dis_data$abund_per_band[data_indices, ]
  bands <- dis_data$bands_per_sample[data_indices]
  max_dist <- dis_data$max_dist[data_indices, ] / 100
  species <- dis_data$species[data_indices]
  max_intervals <- dis_data$max_intervals
  
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

sp_sample_size <- data.frame(table(dis_data$sp_list))
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

(lppd_difference_plot <- ggplot(data = lppd_summary, aes(x = Species, y = Mean_Difference)) +
    geom_errorbar(aes(ymin = Mean_Difference - StdErr, ymax = Mean_Difference + StdErr)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    #coord_flip() +
    NULL)

(lppd_diff_vs_n <- ggplot(data = lppd_summary, aes(x = N, y = Mean_Difference)) +
    geom_point(aes(color = Model_Preference)) +
    stat_smooth(method = "lm", 
                formula = y ~ x, 
                geom = "smooth") +  
    NULL
)

####### Output ####################################

write.table(lppd, file = "data/generated/lppd_distance.csv", sep = ",", row.names = FALSE)
write.table(lppd_summary, file = "data/generated/lppd_summary_distance.csv", sep = ",", row.names = FALSE)

png(filename = "output/plots/cv_distance.png",
    width = 20, height = 6, units = "in", res = 300)
print(lppd_difference_plot)
dev.off()

png(filename = "output/plots/cv_distance_diff_vs_n.png",
    width = 10, height = 6, units = "in", res = 300)
print(lppd_diff_vs_n)
dev.off()
