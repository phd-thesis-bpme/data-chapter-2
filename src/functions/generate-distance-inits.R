generate_distance_inits <- function(n_chains = NULL,
                                    napops_skip = NULL,
                                    dis_data = NULL)
{
  library(napops)
  
  binomial <- read.csv("data/generated/binomial_names.csv")
  
  sp <- data.frame(Scientific_BT = gsub(pattern = "_",
             replacement= " ",
             x = dimnames(dis_data$phylo_corr)[[1]]))
  sp <- dplyr::left_join(sp, binomial, by = "Scientific_BT")
  sp <- sp[, c("Scientific_BT", "Code")]
  sp$edr <- NA
  for (i in 1:nrow(sp))
  {
    if (sp$Code[i] %in% napops_skip)
    {
      sp$edr[i] <- NA
    }else
    {
      sp$edr[i] <- edr(species = sp$Code[i], model = 1, road = TRUE, forest = 1)$EDR_est
    }
  }
  sp[which(is.na(sp$edr)), "edr"] <- mean(sp$edr, na.rm = TRUE)

  inits <- list()
  for (i in 1:n_chains)
  {
    intercept <- rnorm(1, mean = 0.05, sd = 0.1)
    mu_mig_strat <- rnorm(2, mean = 0, sd = 0.05)
    mu_habitat <- rnorm(2, mean = 0, sd = 0.05)
    beta_mass <- rnorm(1, mean = 0.01, sd = 0.005)
    beta_pitch <- rnorm(1, mean = -0.01, sd = 0.005)
    sigma <- rexp(1, rate = 5)
    
    log_tau_mean <- log(sp$edr / 100)
    log_tau <- MASS::mvrnorm(n = 1, mu = log_tau_mean, Sigma = diag(0.001, length(log_tau_mean)))
    log_tau[which(log_tau >= 1.5)] <- 1.4999
    log_tau[which(log_tau <= -2)] <- -1.9999
    
    inits[[i]] <- list(intercept = intercept,
                       mu_mig_strat = mu_mig_strat,
                       mu_habitat = mu_habitat,
                       beta_mass = beta_mass,
                       beta_pitch = beta_pitch,
                       sigma = sigma,
                       log_tau = log_tau)
  }
  
  return(inits)
}
