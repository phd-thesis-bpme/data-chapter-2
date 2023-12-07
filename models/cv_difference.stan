data {
  int<lower=0> N;
  int<lower=0> n_species;
  array [N] real difference; //response variable
  array [N] int species;
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

