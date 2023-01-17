####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 06-removal-analysis.R
# Created April 2022
# Last Updated January 2023

####### Import Libraries and External Files #######

library(cmdstanr)
library(napops)
library(plyr)
library(ggplot2)

####### Read Data #################################

rem_model <- readRDS("output/model_runs/removal_predictions.RDS")
load("data/generated/corr_matrix_predict.rda")
binomial <- read.csv("data/generated/binomial_names.csv")
load("data/generated/removal_stan_data_pred.rda")

####### Data Wrangling ############################

# Extract log_phi summary statistics from full Stan model runs
rem_summary <- rem_model$summary(variables = "log_phi")

# Add species names to these summaries
rem_summary$Scientific_BT <- gsub("_", " ", rownames(corr_matrix_predict))
rem_summary <- join(rem_summary, binomial[, c("Scientific_BT", "Code")], by = "Scientific_BT")

# Get data sample size for all species and add to summary
species_n <- data.frame(table(removal_stan_data_pred$species))
names(species_n) <- c("index", "N")
rem_summary$index <- seq(1, nrow(rem_summary))
rem_summary$N <- 0
for (i in 1:nrow(species_n))
{
  rem_summary[which(rem_summary$index == species_n$index[i]), "N"] <-
    species_n$N[i]
}

# Get original single-species NA-POPS estimates
napops_summary <- napops::coef_removal(species = rem_summary$Code, model = 1)

####### 1-to-1 Plot ###############################

to_plot <- merge(rem_summary[,c("Code", "mean")], napops_summary[, c("Species", "Intercept")],
                 by.x = "Code", by.y = "Species")
names(to_plot) <- c("Species", "Multi", "Single")
to_plot$Label <- ""
to_plot$diff <- NA
for (i in 1:nrow(to_plot))
{
  to_plot$diff[i] <- abs(to_plot$Multi[i] - to_plot$Single[i])
  if (to_plot$diff[i] > 0.1)
  {
    to_plot$Label[i] <- to_plot$Species[i]
  }
}
single_vs_multi_plot <- ggplot(data = to_plot, mapping = aes(x = exp(Single), y = exp(Multi))) +
  geom_point() +
  geom_abline(slope = 1) +
  xlim(0,1) +
  ylim(0,1) +
  geom_text(aes(label = Label)) +
  NULL

####### SD Comparison Plot ########################

# need the napops R package to allow user to retrieve variation

####### Species-specific Plots ####################

sp <- c("LCTH", "LEPC", "HASP", "TRBL", "SPOW", "KIWA", "BITH")

to_plot <- rem_summary[which(rem_summary$Code %in% sp), ]

species_plot <- ggplot(data = to_plot, aes(x = Code, y = exp(mean),
                                           ymin = exp(q5), ymax = exp(q95))) +
  geom_point() + 
  geom_errorbar() +
  xlab("Species") +
  ylab("Cue Rate") +
  geom_text(aes(y = exp(q95) + 0.25, label = N)) +
  NULL

sp_plot_list <- vector(mode = "list", length = length(sp))
names(sp_plot_list) <- sp

for (s in sp)
{
  to_plot <- rem_summary[which(rem_summary$Code == s), ]
  sp_plot_list[[sp]] <- ggplot(to_plot, aes())
}

####### Output ####################################

png("plots/removal_comparison.png",
    width = 6, height = 4, res = 600, units = "in")
print(comparison_plot)
dev.off()
