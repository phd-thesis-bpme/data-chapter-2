data {
  int<lower=0> N;
  int<lower=0> n_species;
  vector[N] difference; //response variable
  int species[N];
}

parameters {
  real overall_difference;
  vector[n_species] species_difference;
  //real<lower = 1> nu;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] mu;
  
  mu = overall_difference + species_difference[species];
}

model {
  //difference ~ student_t(nu, mu, sigma);
  difference ~ normal(mu, sigma);
  
  //nu ~ gamma(2, 0.1);
  overall_difference ~ std_normal();
  species_difference ~ std_normal();
  sigma ~ exponential(1);
}

