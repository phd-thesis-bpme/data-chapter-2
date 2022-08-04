####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 02-removal-model.R
# Created April 2022
# Last Updated August 2022

####### Import Libraries and External Files #######

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

####### Read Data #################################

load("data/generated/removal_stan_data.rda")

####### Stan Modelling ############################

model <- stan_model(file = "models/removal.stan")

removal_stan_job <- sampling(model,
                            data = removal_stan_data,
                            verbose = TRUE,
                            chains = 4,
                            iter = 1000,
                            warmup = 500,
                            cores = 4,
                            pars = c("log_phi", "sigma", "mu"),
                            control = list(adapt_delta = 0.8,
                                           max_treedepth = 15))

####### Output ####################################

save(removal_stan_job, file = "data/generated/removal_model.rda")

