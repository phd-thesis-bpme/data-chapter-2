functions {
  real partial_sum_lpmf(int [, ] slice_abund_per_band,
                        int start, int end,
                        int max_intervals,
                        int [] bands_per_sample,
                        int [, ] max_dist,
                        int [] species,
                        vector log_tau)
  {
    real lp = 0;
    int Pi_size = end - start + 1;
    int Pi_index = 1;
    matrix[Pi_size, max_intervals] Pi = rep_matrix(0, Pi_size, max_intervals);
    
    for (i in start:end)
    {
      for (k in 2:bands_per_sample[i])
      {
        Pi[Pi_index,k] = ((1 - exp(-(max_dist[i,k]^2 / exp(log_tau[species[i]])^2))) - 
        (1 - exp(-(max_dist[i,k - 1]^2 / exp(log_tau[species[i]])^2)))) / 
        (1 - exp(-(max_dist[i,bands_per_sample[i]]^2 / exp(log_tau[species[i]])^2)));
      }
      Pi[Pi_index,1] = 1 - sum(Pi[Pi_index,]);
      
      lp = lp + multinomial_lpmf(slice_abund_per_band[Pi_index, ] | to_vector(Pi[Pi_index, ]));
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
  int<lower = 1> n_species;           // total number of specise
  int<lower = 2> max_intervals;       // maximum number of intervals being considered
  int<lower = 1> grainsize;           // grainsize for reduce_sum() function
  int species[n_samples];             // species being considered for each sample
  int abund_per_band[n_samples, max_intervals];// abundance in distance band k for sample i
  int bands_per_sample[n_samples]; // number of distance bands for sample i
  int max_dist[n_samples, max_intervals]; // max distance for distance band k
  corr_matrix[n_species] phylo_corr; // correlation matrix of phylogeny
  real<lower = 0> lambda;            //Pagel's lambda
}

transformed data {
  corr_matrix[n_species] phylo_corr_pl;
  phylo_corr_pl = phylo_pl(phylo_corr, lambda);
  
}

parameters {
  row_vector<lower = 3>[n_species] mu;
  real<lower = 0> sigma;
  vector<lower = 0>[n_species] log_tau;
}

model {
  sigma ~ exponential(1);
  mu ~ lognormal(2, 0.5);

 // mu_vector = rep_array(mu, n_species);
  
  log_tau ~ multi_normal(mu, quad_form_diag(phylo_corr_pl, rep_vector(sigma, n_species)));

  target += reduce_sum(partial_sum_lpmf,
                       abund_per_band,
                       grainsize,
                       max_intervals,
                       bands_per_sample,
                       max_dist,
                       species,
                       log_tau);
}
