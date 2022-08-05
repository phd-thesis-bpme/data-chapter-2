####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 02-removal-model.R
# Created April 2022
# Last Updated August 2022

####### Import Libraries and External Files #######

library(cmdstanr)
#options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)

####### Read Data #################################

load("data/generated/removal_stan_data.rda")

####### Stan Modelling ############################

model <- cmdstan_model(stan_file = "models/removal.stan",
                       cpp_options = list(stan_threads = TRUE))#model <- stan_model(file = "models/removal.stan")

removal_stan_data$grainsize <- 1

removal_stan_fit <- model$sample(
  data = removal_stan_data,
  chains = 4,
  parallel_chains = 4,
  refresh = 100,
  threads_per_chain = 4
)

# removal_stan_job <- sampling(model,
#                             data = removal_stan_data,
#                             verbose = TRUE,
#                             chains = 4,
#                             iter = 1000,
#                             warmup = 500,
#                             cores = 4,
#                             pars = c("log_phi", "sigma", "mu"),
#                             control = list(adapt_delta = 0.8,
#                                            max_treedepth = 15))

####### Output ####################################

removal_stan_fit$save_object(file = "data/generated/removal_model.RDS")
#save(removal_stan_fit, file = "data/generated/removal_model.rda")

