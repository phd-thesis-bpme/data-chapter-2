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

####### Data Wrangling ############################

# Extract log_phi summary statistics from full Stan model runs
rem_summary <- rem_model$summary(variables = "log_phi")

# Add species names to these summaries
rem_summary$Scientific_BT <- gsub("_", " ", rownames(corr_matrix_predict))
rem_summary <- join(rem_summary, binomial[, c("Scientific_BT", "Code")], by = "Scientific_BT")

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


####### Output ####################################

png("plots/removal_comparison.png",
    width = 6, height = 4, res = 600, units = "in")
print(comparison_plot)
dev.off()
