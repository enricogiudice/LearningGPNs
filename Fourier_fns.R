# These are the functions to generate non-linear data with the Fourier transform approach
# lambda=0 is linear, larger is more non-linear

Fou_trans <- function(x, lambda, concentration, n_Four = 6){
  
  n_Four <- min(n_Four, 10) # check we don't have too many
  n_Four <- max(n_Four, 2) # or too few!
  
  concentration <- max(concentration, 1e-3) # check it remains positive
  
  # store the multiples of x
  x_mults <- matrix(rep(x, each = n_Four)*(1:n_Four), nrow = n_Four)
  
  # compute the weights for the frequencies (we have most twice)
  if(lambda == 0) {
    freq_weights <- c(1, rep(0, n_Four, each = 2))
  }
  
  else {
    freq_weights <- exp(-rep(0:n_Four, each = 2)[-1] / lambda)
    freq_weights <- freq_weights/sum(freq_weights)
  }
  
  # Dirichlet sampler
  sampled_weights <- rgamma(length(freq_weights), freq_weights*concentration)
  sampled_weights <- sampled_weights/sum(sampled_weights)
  # chose random signs
  sampled_weights <- sampled_weights*sample(c(-1, 1), length(sampled_weights), replace = TRUE)
  # add beta
  sampled_weights <- runif(1, 0.5, 2)*sampled_weights
  
  # transform x by taking cosine terms, weighted and ignoring intercept
  x_transform_cos <- colSums(sampled_weights[2*1:n_Four]*cos(x_mults))
  # transform x by taking sine terms, weighted and ignoring intercept
  x_transform_sin <- colSums(sampled_weights[2*1:n_Four + 1]*sin(x_mults))  
  # include linear term  
  x_transform <- sampled_weights[1]*x + x_transform_cos + x_transform_sin
  
  return(x_transform) 
}

Fou_nldata <- function(dag, N, lambda, noise.sd, concentration = 1, standardize = FALSE) {
  n <- ncol(dag)  # number of variables
  data <- matrix(0, nrow = N, ncol = n)  # to store the simulated data
  top_order <- rev(BiDAG:::DAGtopartition(n, dag)$permy)  # go down order
  
  for (jj in top_order) {
    parents <- which(dag[, jj] == 1)  # find parents
    lp <- length(parents)  # number of parents
    
    if (lp == 0) {  # zero parents
      data[ ,jj] <- rnorm(N, 0, noise.sd)
    }
    else if (lp == 1) {  # one parent
      trans.pa <- Fou_trans(data[ ,parents], lambda, concentration)
      data[ ,jj] <- trans.pa + rnorm(N, 0, noise.sd)
    }
    else {  # More than one parent
      trans.pa <- Fou_trans(data[ ,parents], lambda, concentration)
      data[, jj] <- rowSums(trans.pa) + rnorm(N, 0, noise.sd)
    }
  }
  
  if(standardize) { 
    data <- scale(data)
  } 
  return(data)
}
