data {
  int<lower = 1> n_samples;           // total number of sampling events i
  int<lower = 1> n_covariates;        // total number of covariates
  int<lower = 2> max_intervals;       // maximum number of intervals being considered
  int abund_per_band[n_samples, max_intervals];// abundance in distance band k for sample i
  vector[n_samples] abund_per_sample; // total abundnace for sample i
  int bands_per_sample[n_samples]; // number of distance bands for sample i
  matrix[n_samples, max_intervals] max_dist; // max distance for distance band k
  matrix[n_samples, n_covariates] X;  // matrix of covariates
}

parameters {
  vector<lower=1>[n_covariates] theta;
}

model {
  vector[n_samples] log_tau;             // effective detection radius
  matrix[n_samples, max_intervals] Pi;   // probabilities
  
  theta ~ lognormal(1, 1);
  log_tau = X * theta;
  
  Pi = rep_matrix(0, n_samples, max_intervals);
  
  for (i in 1:n_samples)
  {
    for (k in 2:bands_per_sample[i])
    {
      Pi[i,k] = ((1 - exp(-(max_dist[i,k]^2 / exp(log_tau[i])^2))) - 
      (1 - exp(-(max_dist[i,k - 1]^2 / exp(log_tau[i])^2)))) / 
      (1 - exp(-(max_dist[i,bands_per_sample[i]]^2 / exp(log_tau[i])^2)));
    }
    Pi[i,1] = 1 - sum(Pi[i,]);
    
    abund_per_band[i,] ~ multinomial(to_vector(Pi[i,]));
  }
}
