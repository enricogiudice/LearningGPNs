data {
  int<lower=1> N_obs;
  vector[N_obs] y_obs;
}

parameters {
  real<lower=-2, upper=2> mu;
  real<lower=0, upper=2> sigma;
}

model {
  target += std_normal_lpdf(mu);
  target += inv_gamma_lpdf(sigma | 1, 1);
  target += normal_lpdf(y_obs | mu, sigma);
}
