data {
  int<lower=1> N_obs;
  int<lower=1> d;
  real X[N_obs, d];
  vector[N_obs] y_obs;
}

parameters {
  vector<lower=0, upper=10>[d] rho;
  real<lower=0> var1;
  real<lower=0> var2;
  real<lower=-2, upper=2> mu;
  real<lower=0, upper=2> sigma;
}

model {
  matrix[N_obs, N_obs] e1 = rep_matrix(0, N_obs, N_obs);
  matrix[N_obs, N_obs] e2 = rep_matrix(0, N_obs, N_obs);
  matrix[N_obs, N_obs] s2 = rep_matrix(0, N_obs, N_obs);
  matrix[N_obs, N_obs] cov;
  matrix[N_obs, N_obs] L_cov;
  
  for (i in 1:d) {
    cov = cov_exp_quad(X[ ,i], var1, rho[i]);
    e1 += cov;
    s2 += cov .* cov;
  }
  e2 = 0.5 * var2 * (e1 .* e1 - s2) + e1;
  
  cov = e2 + diag_matrix(rep_vector(square(sigma), N_obs));
  L_cov = cholesky_decompose(cov);
  
  target += inv_gamma_lpdf(rho | 2, 2);
  target += inv_gamma_lpdf(var1 | 1, 1);
  target += inv_gamma_lpdf(var2 | 1, 1);
  target += std_normal_lpdf(mu);
  target += inv_gamma_lpdf(sigma | 1, 1);
  target += multi_normal_cholesky_lpdf(y_obs | rep_vector(mu, N_obs), L_cov);
}
