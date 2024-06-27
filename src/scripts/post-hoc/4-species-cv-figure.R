####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# posthoc/4-species-cv-figure.R
# Created April 2024
# Last Updated June 2024

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

load("data/generated/distance_stan_data_cv.rda")
distance_all <- readRDS("output/model_runs/distance_ms.RDS")
distance_folds <- vector(mode = "list", length = k)
for (i in 1:k)
{
  distance_folds[[i]] <- readRDS(paste0("output/model_runs/cv_distance/species/fold_", i, ".RDS"))
}
dis_fold_nums <- read.csv("data/generated/distance_cv_folds_species.csv")

traits <- read.csv("data/raw/traits.csv")

####### Removal Model #################################

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
rem_cv_df$Percent_Diff <- rem_cv_df$Difference / rem_cv_df$True_Exp

rem_model_data <- list(n_samples = nrow(rem_cv_df),
                       true_values = rem_cv_df$True_Exp,
                       predicted_values = rem_cv_df$Predicted_Exp)

rem_cv_model <- cmdstan_model(stan_file = "models/predicted-vs-true.stan")
rem_cv_model_run <- rem_cv_model$sample(data = rem_model_data,
                                          iter_warmup = 1000,
                                          iter_sampling = 2000,
                                          chains = 4,
                                          parallel_chains = 4,
                                          refresh = 10)
rem_cv_model_draws <- rem_cv_model_run$draws(format = "df")
rem_cv_model_summary <- rem_cv_model_run$summary()

removal_cv_plot <- ggplot(data = rem_cv_df, aes(x = True_Exp, y = Predicted_Exp)) +
  geom_abline(intercept = rem_cv_model_draws$intercept, slope = rem_cv_model_draws$slope,
              color = "grey", alpha = 0.1) +
  geom_abline(intercept = mean(rem_cv_model_draws$intercept), slope = mean(rem_cv_model_draws$slope),
              color = "black") +
  geom_point() +
  geom_abline(slope = 1, color = "red", linetype = 2) +
  xlim(0,0.8) + ylim(0,0.8) +
  xlab("Cue Rate (Full Data)") +
  ylab("Cue Rate (No Data)") +
  NULL

####### Distance Model #################################

tau_all <- distance_all$summary("log_tau")
tau_all$Species <- unique(distance_stan_data_cv$sp_list)

dis_cv_df <- dis_fold_nums[!duplicated(dis_fold_nums$Species), ]
dis_cv_df <- merge(dis_cv_df, tau_all[, c("mean", "Species")],
                   by = "Species")
dis_cv_df$Predicted <- NA
names(dis_cv_df)[3] <- "True"

for (i in 1:k)
{
  tau_fold <- distance_folds[[i]]$summary("log_tau")
  tau_fold$Species <- unique(distance_stan_data_cv$sp_list)
  
  sps <- dis_cv_df[which(dis_cv_df$cv_fold == i), "Species"]
  
  for (sp in sps)
  {
    dis_cv_df[which(dis_cv_df$Species == sp), "Predicted"] <- tau_fold[which(tau_fold$Species == sp), "mean"]
  }
}

dis_cv_df$True_Exp <- exp(dis_cv_df$True) * 100
dis_cv_df$Predicted_Exp <- exp(dis_cv_df$Predicted) * 100
dis_cv_df$Difference <- dis_cv_df$True_Exp - dis_cv_df$Predicted_Exp
dis_cv_df$Abs_Diff <- abs(dis_cv_df$Difference)

dis_model_data <- list(n_samples = nrow(dis_cv_df),
                       true_values = dis_cv_df$True_Exp,
                       predicted_values = dis_cv_df$Predicted_Exp)

dis_cv_model <- cmdstan_model(stan_file = "models/predicted-vs-true.stan")
dis_cv_model_run <- dis_cv_model$sample(data = dis_model_data,
                                        iter_warmup = 1000,
                                        iter_sampling = 2000,
                                        chains = 4,
                                        parallel_chains = 4,
                                        refresh = 10)
dis_cv_model_draws <- dis_cv_model_run$draws(format = "df")
dis_cv_model_summary <- dis_cv_model_run$summary()

distance_cv_plot <- ggplot(data = dis_cv_df, aes(x = True_Exp, y = Predicted_Exp)) +
  geom_abline(intercept = dis_cv_model_draws$intercept, slope = dis_cv_model_draws$slope,
              color = "grey", alpha = 0.1) +
  geom_abline(intercept = mean(dis_cv_model_draws$intercept), slope = mean(dis_cv_model_draws$slope),
              color = "black") +
  geom_point() +
  geom_abline(slope = 1, color = "red", linetype = 2) +
  xlim(0,500) + ylim(0,500) +
  xlab("EDR (Full Data)") +
  ylab("EDR (No Data)") +
  NULL

dis_cv_df <- merge(dis_cv_df, traits[, c("Code", "Migrant", "Habitat")],
                   by.x = "Species", by.y = "Code")
dis_cv_df$Migrant <- ifelse(dis_cv_df$Migrant == 1,
                            "Migrant",
                            "Resident")
dis_cv_df$Habitat <- ifelse(dis_cv_df$Habitat == 1,
                            "Closed",
                            "Open")
dis_cv_df$Trait_Group <- paste0(dis_cv_df$Habitat, "-", dis_cv_df$Migrant)

distance_cv_traits_plot <- ggplot(data = dis_cv_df, aes(x = True_Exp, y = Predicted_Exp)) +
  geom_abline(intercept = dis_cv_model_draws$intercept, slope = dis_cv_model_draws$slope,
              color = "grey", alpha = 0.1) +
  geom_abline(intercept = mean(dis_cv_model_draws$intercept), slope = mean(dis_cv_model_draws$slope),
              color = "black") +
  geom_point(aes(color = Trait_Group)) +
  geom_abline(slope = 1, color = "red", linetype = 2) +
  xlim(60,150) + ylim(60,150) +
  theme(legend.position = "right") +
  NULL


####### Output ####################################

write.table(rem_cv_model_summary, file = "data/generated/rem_species_cv_model.csv",
            sep = ",", row.names = FALSE)
write.table(dis_cv_model_summary, file = "data/generated/dis_species_cv_model.csv",
            sep = ",", row.names = FALSE)

tiff("output/plots/species_cv_plot.tiff",
    width = 6, height = 3, res = 300, units = "in")
ggarrange(removal_cv_plot, distance_cv_plot, ncol = 2, labels = c("A", "B"))
dev.off()

png("output/plots/species_cv_plot.png",
    width = 6, height = 3, res = 300, units = "in")
ggarrange(removal_cv_plot, distance_cv_plot, ncol = 2, labels = c("A", "B"))
dev.off()

tiff("output/plots/species_cv_distance_traits.tiff",
     width = 5, height = 3, res = 300, units = "in")
print(distance_cv_traits_plot)
dev.off()

png("output/plots/species_cv_distance_traits.png",
    width = 5, height = 3, res = 300, units = "in")
print(distance_cv_traits_plot)
dev.off()
