####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 7-removal-cv-analysis.R
# Created October 2023
# Last Updated October 2023

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

####### Read Data #################################

load("data/generated/removal_stan_data_cv.rda")
cv_folds <- read.csv("data/generated/removal_cv_folds.csv")
phylo_tree <- ape::read.nexus(file = "data/raw/all_species.nex")
binomial <- read.csv("data/generated/binomial_names.csv")
traits <- read.csv("data/raw/traits.csv")

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

lppd_summary <- merge(lppd_summary, binomial[,c("Code", "Scientific_BT")],
                      by.x = "Species", by.y = "Code")
lppd_summary$Scientific_BT <- gsub(" ", "_", lppd_summary)
lppd_summary <- merge(lppd_summary, traits[, c("Code", "Migrant")],
                      by.x = "Species", by.y = "Code")
tree <- ape::consensus(phylo_tree)
drops <- binomial[-which(binomial$Code %in% removal_stan_data_cv$sp_list$Species), "Scientific_BT"]
drops <- gsub(" ", "_", drops)
tree <- drop.tip(tree, drops)


(model_pref_tree <- ggtree(tree, branch.length="none") %<+% lppd_summary[, c("Scientific_BT", "Model_Preference", "N", "Migrant")] +
    geom_tiplab(offset = .6, hjust = .5) +
  geom_tippoint(aes(color = Model_Preference, shape = as.factor(Migrant), size = log(N))))


# cv_model <- cmdstan_model(stan_file = "models/cv_difference.stan")
# cv_model_data <- list(N = nrow(lppd),
#                       n_species = length(unique(lppd$Species)),
#                       difference = lppd$Difference,
#                       species = removal_stan_data_cv$species)
# 
# cv_model_run <- cv_model$sample(
#   data = cv_model_data,
#   iter_warmup = 25,
#   iter_sampling = 25,
#   chains = 1,
#   parallel_chains = 4,
#   refresh = 10
# )
# 
# saveRDS(cv_model_run, file = "output/model_runs/cv_removal/cv_model_removal.RDS")
# 
# cv_model_summary <- cv_model_run$summary(c("nu", "overall_difference", "species_difference", "sigma"))

####### Output ####################################

write.table(lppd, file = "data/generated/lppd_removal.csv", sep = ",", row.names = FALSE)
write.table(lppd_summary, file = "data/generated/lppd_summary_removal.csv", sep = ",", row.names = FALSE)

png(filename = "output/plots/cv_removal.png",
    width = 20, height = 6, units = "in", res = 300)
print(lppd_difference_plot)
dev.off()

png(filename = "output/plots/cv_removal_diff_vs_n.png",
    width = 10, height = 6, units = "in", res = 300)
print(lppd_diff_vs_n)
dev.off()

png(filename = "output/plots/cv_removal_tree.png",
    width = 30, height = 60, units = "in", res = 300)
print(model_pref_tree)
dev.off()

