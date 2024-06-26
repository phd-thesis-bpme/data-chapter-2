
data {
  int<lower = 0> n_samples;
  vector [n_samples] true_values;
  vector [n_samples] predicted_values;
}

parameters {
  real intercept;
  real slope;
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  predicted_values ~ normal(intercept + slope * true_values, sigma);
  intercept ~ normal(0,1);
  slope ~ normal(1,1);
  sigma ~ exponential(1);
}

