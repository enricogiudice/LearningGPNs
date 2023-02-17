library(rstan)
library(bridgesampling)
options(mc.cores = parallel::detectCores())

# Compile stan files
Gauss_mod <- stan_model("/Users/giudic0000/Downloads/Nonlinear scoring/Structure Learning/Gauss.stan") 
GP_mod <- stan_model("/Users/giudic0000/Downloads/Nonlinear scoring/Structure Learning/Add.stan") 

# Laplace approximation of log marginal likelihood
Laplace <- function(optim) {
  d = length(optim$par)
  penalty <- determinant(-optim$hessian, logarithm = T)$modulus - d*log(2*pi)
  optim$value - 0.5*penalty
}

# Optimizes posterior of Gaussian and computes Laplace approximation
Gauss.lap <- function(y) {  
  Sdata <- list(N_obs = length(y), y_obs = y)
  init <- list(mu = 0, sigma = 1)
  
  opt <- optimizing(Gauss_mod, data = Sdata, init = init, hessian = T)
  Laplace(opt)
}

# Optimizes posterior of additive GP
GP.lap <- function(y, X) {  
  X <- as.matrix(X)
  d <- ncol(X)
  Sdata <- list(N_obs = nrow(X), d = d, X = X, y_obs = y)  
  init <- list(rho = as.array(rep(1, d)), mu = 0, sigma = 0.5)
  
  opt <- optimizing(GP_mod, data = Sdata, init = init, hessian = T)
  Laplace(opt)
}

# Bridgesampling for marginal likelihoods
Gauss.mcmc <- function(y) {
  N <- length(y)
  stanfit <- sampling(Gauss_mod, data = list(N_obs = N, y_obs = y), 
                      iter = 500, warmup = 200, chains = 1, refresh = 0)
  bridge_sampler(stanfit, silent = TRUE)$logml
}

GP.mcmc <- function(y, X) { 
  N <- length(y)
  X <- as.matrix(X)
  d <- ncol(X)
  stanfit <- sampling(GP_mod, data = list(N_obs = N, d = d, X = X, y_obs = y),  
                      iter = 500, warmup = 200, chains = 1, refresh = 0)
  bridge_sampler(stanfit, silent = TRUE)$logml
}

