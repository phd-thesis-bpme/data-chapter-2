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
  int<lower = 1> n_mig_strat;        //total number of migration strategies
  int mig_strat[n_species];        //migration strategy for each species
  int<lower = 1> n_habitat;        //total number of habitat preferences
  int habitat[n_species];        //habitat preference for each species
  real<lower = 0> mass[n_species];  //log mass of species
  real<lower = 0> pitch[n_species]; //song pitch of species
}

transformed data {
  corr_matrix[n_species] phylo_corr_pl;
  phylo_corr_pl = phylo_pl(phylo_corr, lambda);
  
}

parameters {
  row_vector<lower = 3>[n_species] mu;
  row_vector[n_mig_strat] mu_mig_strat;
  row_vector[n_habitat] mu_habitat;
  real beta_mass;
  real beta_pitch;
  real<lower = 0> sigma;
  vector<lower = 0>[n_species] log_tau;
}

model {
  log_tau ~ multi_normal(mu, quad_form_diag(phylo_corr_pl, rep_vector(sigma, n_species)));
  
  for (sp in 1:n_species)
  {
    mu[sp] ~ lognormal(mu_mig_strat[mig_strat[sp]] +
                       mu_habitat[habitat[sp]] +
                       beta_mass * mass[sp] +
                       beta_pitch * pitch[sp], 
                       0.5);
  }
  mu_mig_strat ~ normal(0,1);
  mu_habitat ~ normal(0,1);
  beta_mass ~ normal(0,1);
  beta_pitch ~ normal(0,1);
  
  sigma ~ exponential(1);

  target += reduce_sum(partial_sum_lpmf,
                       abund_per_band,
                       grainsize,
                       max_intervals,
                       bands_per_sample,
                       max_dist,
                       species,
                       log_tau);
}
