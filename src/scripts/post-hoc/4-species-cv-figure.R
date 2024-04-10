####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# posthoc/4-species-cv-figure.R
# Created April 2024
# Last Updated April 2024

####### Import Libraries and External Files #######

library(cmdstanr)
library(ggplot2)
library(ggpubr)
theme_set(theme_pubclean())

####### Set Constants #############################

k <- 10

####### Read Data #################################

load("data/generated/removal_stan_data_cv.rda")
removal_all <- readRDS("output/model_runs/removal_ms.RDS")
removal_folds <- vector(mode = "list", length = k)
for (i in 1:k)
{
  removal_folds[[i]] <- readRDS(paste0("output/model_runs/cv_removal/species/fold_", i, ".RDS"))
}
rem_fold_nums <- read.csv("data/generated/removal_cv_folds_species.csv")

####### Main Code #################################

phi_all <- removal_all$summary("log_phi")
phi_all$Species <- removal_stan_data_cv$sp_all

rem_cv_df <- rem_fold_nums[!duplicated(rem_fold_nums$Species), ]
rem_cv_df <- merge(rem_cv_df, phi_all[, c("mean", "Species")],
                   by = "Species")
rem_cv_df$Predicted <- NA
names(rem_cv_df)[3] <- "True"

for (i in 1:k)
{
  phi_fold <- removal_folds[[i]]$summary("log_phi")
  phi_fold$Species <- removal_stan_data_cv$sp_all
  
  sps <- rem_cv_df[which(rem_cv_df$cv_fold == i), "Species"]
  
  for (sp in sps)
  {
    rem_cv_df[which(rem_cv_df$Species == sp), "Predicted"] <- phi_fold[which(phi_fold$Species == sp), "mean"]
  }
}

rem_cv_df$True_Exp <- exp(rem_cv_df$True)
rem_cv_df$Predicted_Exp <- exp(rem_cv_df$Predicted)
rem_cv_df$Difference <- rem_cv_df$True_Exp - rem_cv_df$Predicted_Exp
rem_cv_df$Abs_Diff <- abs(rem_cv_df$Difference)

(removal_cv_plot <- ggplot(data = rem_cv_df, aes(x = True_Exp, y = Predicted_Exp)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0,0.8) + ylim(0,0.8) +
  NULL)

####### Output ####################################