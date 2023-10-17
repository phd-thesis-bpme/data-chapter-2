data {
  int<lower=0> N;
  int<lower=0> n_species;
  vector[N] difference; //response variable
  int species[N];
}

parameters {
  real overall_difference;
  vector[n_species] species_difference;
  real<lower=0> sigma;
}

model {
  difference ~ normal(overall_difference + species_difference[species], sigma);
  
  overall_difference ~ normal(0,0.1);
  species_difference ~ normal(0,0.1);
  sigma ~ exponential(5);
}

