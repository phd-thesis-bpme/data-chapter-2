functions {
  real partial_sum_lpmf(int [, ] slice_abund_per_band,
                        int start, int end,
                        int max_intervals,
                        int [] bands_per_sample,
                        int [, ] max_time,
                        int [] species,
                        vector log_phi)
  {
    real lp = 0;
    int Pi_size = end - start + 1;
    int Pi_index = 1;
    matrix[Pi_size, max_intervals] Pi = rep_matrix(0, Pi_size, max_intervals);   // probabilities

    for (i in start:end)
    {
      for (j in 2:bands_per_sample[i])
      {
        Pi[Pi_index,j] = (exp(-max_time[i,j-1] * exp(log_phi[species[i]])) -
                   exp(-max_time[i,j] * exp(log_phi[species[i]]))) /
                  (1 - exp(-max_time[i,bands_per_sample[i]] * exp(log_phi[species[i]])));
      }
      Pi[Pi_index,1] = 1 - sum(Pi[Pi_index,]);
      
      lp = lp + multinomial_lpmf(slice_abund_per_band[Pi_index,] | to_vector(Pi[Pi_index,]));
      Pi_index = Pi_index + 1;
    }

    return lp;
  }
  
   // Adapted from Statistical Rethinking Ed 2, Chapter 14 Page 484 
  matrix phylo_pl(matrix x,
                  real lambda)
  {
    int N = dims(x)[1];
    matrix[N,N] K;
    
    K = x * lambda;
    
    for (i in 1:N)
    {
      K[i,i] = 1;
    }
    
    return K;
  }

}

data {
  int<lower = 1> n_samples;           // total number of sampling events i
  int<lower = 2> max_intervals;       // maximum number of intervals being considered
  int<lower = 1> n_species;           // total number of species being modelled
  int<lower = 1> grainsize;           // grainsize for reduce_sum() function
  int species[n_samples];             // species being considered for each sample
  int abund_per_band[n_samples, max_intervals];// abundance in time band j for sample i
  int bands_per_sample[n_samples]; // number of time bands for sample i
  int max_time[n_samples, max_intervals]; // max time duration for time band j
  corr_matrix[n_species] phylo_corr; // correlation matrix of phylogeny
  real<lower = 0> lambda;            //Pagel's lambda
  int<lower = 1> n_mig_strat;        //total number of migration strategies
  int mig_strat[n_species];        //migration strategy for each species
}

transformed data {
  corr_matrix[n_species] phylo_corr_pl;
  phylo_corr_pl = phylo_pl(phylo_corr, lambda);
}

parameters {
  row_vector[n_species] mu;
  row_vector[n_mig_strat] mu_mig_strat_raw;
  real<lower = 0> sigma;
  vector[n_species] log_phi;
}

transformed parameters {
  row_vector[n_mig_strat] mu_mig_strat;
  mu_mig_strat = -1 + 0.1*mu_mig_strat_raw; //we expect log_phi to be negative
}

model {
  log_phi ~ multi_normal(mu, phylo_corr_pl * sigma);// quad_form_diag(phylo_corr_pl, rep_vector(sigma, n_species)));
  
  for (sp in 1:n_species)
  {
    mu[sp] ~ normal(mu_mig_strat[mig_strat[sp]], 1);
  }
  mu_mig_strat_raw ~ std_normal();
  
  sigma ~ exponential(1);
  
  target += reduce_sum(partial_sum_lpmf,
                       abund_per_band,
                       grainsize,
                       max_intervals,
                       bands_per_sample,
                       max_time,
                       species,
                       log_phi);

}
