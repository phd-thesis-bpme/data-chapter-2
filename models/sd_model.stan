data {
  int<lower = 0> n_samples;
  int<lower = 0> N[n_samples];
  real<lower = 0> stdev[n_samples];
  int<lower = 0> model_type[n_samples];
}

transformed data {
  real log_N[n_samples];
  log_N = log(N);
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real alpha[2];
  real beta[2];

  real<lower=0> sigma;
}

transformed parameters {
  real mu[n_samples];
  for (i in 1:n_samples)
  {
     mu[i] = alpha[model_type[i]] + beta[model_type[i]] * log_N[i]; 
  }

}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  sigma ~ exponential(1);
  alpha ~ normal(0,1);
  beta ~ normal(0,1);
  for (i in 1:n_samples)
  {
    target += normal_lpdf(stdev[i] | mu[i], sigma);
  }

}

