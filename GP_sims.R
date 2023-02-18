library(BiDAG)
library(matrixStats)

source("Fourier_fns.R")
source("BayesStanFns.R")
source("sampling_fns.R")
source("comparison_algs.R")
source("dualPC.R")
insertSource("GPscore.R", package = "BiDAG")

init.seed <- 100
iter <- 100  # number of simulations
lambdas <- c(0, 0.5, 1)  # non-linearity; zero is linear
dual <- T    # use dualPC
n <- 10      # number of nodes
N <- 100     # number of samples
results <- data.frame()

# Parameters for ROC curves
bge.mus <- c(0.01, 0.1, 0.5, 2, 5)
gp.pars <- c(0.1, 0.5, 1, 10, 100) 
dib.pars <- c(1, 1.25, 1.5, 1.75, 2)
pc.pars <- c(0.004, 0.01, 0.05, 0.12, 0.3)

for(lambda in lambdas) {

  for(k in 1:5) {
    bge.par <- bge.mus[k]
    gp.par <- gp.pars[k]
    dib.par <- dib.pars[k]
    pc.par <- pc.pars[k]
    
    for (i in 1:iter) {
      set.seed(init.seed+i)
      
      # Generate DAG & data
      myDAG <- pcalg::randomDAG(n, prob = 0.2, lB = 1, uB = 2) 
      trueDAG <- as(myDAG, "matrix")
      truegraph <- 1*(trueDAG != 0)
      data <- Fou_nldata(truegraph, N, lambda = lambda, noise.sd = 1, standardize = T) 
      
      # Set initial search spaces
      GP.searchspace = set.searchspace(data, dual, "GP", gp.par)
      bge.searchspace = set.searchspace(data, dual, "bge", bge.par)
      
      # Bge score, partition
      bge.fit <- bge.partition.mcmc(bge.searchspace, order = F)
      bge.comp <- compareGraphs(bge.fit, truegraph)
      results <- rbind(results, c(as.character(bge.comp), bge.par, "BGe, partition", lambda))
      
      # Bge score, order
      bge.fit <- bge.partition.mcmc(bge.searchspace, order = T)
      bge.comp <- compareGraphs(bge.fit, truegraph)
      results <- rbind(results, c(as.character(bge.comp), bge.par, "BGe, order", lambda))
      
      # GP score, partition
      gp.fit <- GP.partition.mcmc(data, GP.searchspace, order = F)
      gp.comp <- compareGraphs(gp.fit, truegraph)
      results <- rbind(results, c(as.character(gp.comp), gp.par, "GP, partition", lambda))
      
      # GP score, order
      gp.fit <- GP.partition.mcmc(data, GP.searchspace, order = T)
      gp.comp <- compareGraphs(gp.fit, truegraph)
      results <- rbind(results, c(as.character(gp.comp), gp.par, "GP, order", lambda))
      
      # DIBS+
      dib.fit <- DiBS(data, 1, p = dib.par)
      dib.comp <- compareGraphs(dib.fit, truegraph)
      results <- rbind(results, c(as.character(dib.comp), dib.par, "DiBS+", lambda))
      
      # k-PC, distance correlation
      kPC.dcc.fit <- kPC.dcc.boot(data, nboots = 10, alpha = pc.par)
      kPC.dcc.comp <- compareGraphs(kPC.dcc.fit, truegraph)
      results <- rbind(results, c(as.character(kPC.dcc.comp), pc.par, "kPC-DC", lambda))
      
      # k-PC, HSIC
      kPC.hsic.fit <- kPC.hsic.boot(data, nboots = 10, alpha = pc.par)
      kPC.hsic.comp <- compareGraphs(kPC.hsic.fit, truegraph)
      results <- rbind(results, c(as.character(kPC.hsic.comp), pc.par, "kPC-HSIC", lambda))
    }
  }
}
colnames(results) <- c("ESHD", "TPR", "FPRn", "parameter", "Scorefn", "Lambda")
saveRDS(results, "Results/Sims_Results.rds")
