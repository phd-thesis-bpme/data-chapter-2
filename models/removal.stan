data {
  int<lower = 1> n_samples;           // total number of sampling events i
  int<lower = 2> max_intervals;       // maximum number of intervals being considered
  int<lower = 1> n_species;           // total number of species being modelled
  int species[n_samples];             // species being considered for each sample
  int abund_per_band[n_samples, max_intervals];// abundance in time band j for sample i
  vector[n_samples] abund_per_sample; // total abundnace for sample i
  int bands_per_sample[n_samples]; // number of time bands for sample i
  matrix[n_samples, max_intervals] max_time; // max time duration for time band j
  corr_matrix[n_species] phylo_corr; // correlation matrix of phylogeny
}

parameters {
  row_vector[n_species] mu;
  vector<lower = 0>[n_species] tau;
  vector[n_species] log_phi;
}

model {
  matrix[n_samples, max_intervals] Pi;   // probabilities
  
  tau ~ cauchy(0, 2.5);
  mu ~ normal(0, 1);
  
  log_phi ~ multi_normal(mu, quad_form_diag(phylo_corr, tau));
  
  Pi = rep_matrix(0, n_samples, max_intervals);
  
  for (i in 1:n_samples)
  {
    for (j in 2:bands_per_sample[i])
    {
      Pi[i,j] = (exp(-max_time[i,j-1] * exp(log_phi[species[i]])) - 
                 exp(-max_time[i,j] * exp(log_phi[species[i]]))) / 
                (1 - exp(-max_time[i,bands_per_sample[i]] * exp(log_phi[species[i]])));
    }
    Pi[i,1] = 1 - sum(Pi[i,]);
    
    abund_per_band[i,] ~ multinomial(to_vector(Pi[i,]));
  }

}
