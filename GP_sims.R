library(BiDAG)
library(matrixStats)

source("/Users/giudic0000/Downloads/Nonlinear scoring/Structure Learning/Fourier_fns.R")
source("/Users/giudic0000/Downloads/Nonlinear scoring/Structure Learning/BayesStanFns.R")
source("/Users/giudic0000/Downloads/Nonlinear scoring/Structure Learning/sampling_fns.R")
source("/Users/giudic0000/Downloads/Nonlinear scoring/Structure Learning/comparison_algs.R")
source("/Users/giudic0000/Downloads/Nonlinear scoring/Structure Learning/dualPC.R")
insertSource("~/Downloads/Nonlinear scoring/Structure Learning/GPscore.R", package = "BiDAG")

init.seed <- 100
iter <- 100  # number of simulations
lambdas <- c(0, 0.5, 1)  # non-linearity; zero is linear
dual <- T    # use dualPC
n <- 10      # number of nodes
N <- 100     # number of samples
results <- data.frame()

# Parameters for ROC curves
bge.mus <- c(0.01, 0.1, 0.5, 2, 5)
gp.pars <- c(0.1, 0.5, 1, 10, 100)  # 0.1 is very slow
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
      results <- compare_results(bge.fit, c(bge.par, "BGe, partition", lambda), results, truegraph)
      
      # Bge score, order
      bge.fit <- bge.partition.mcmc(bge.searchspace, order = T)
      results <- compare_results(bge.fit, c(bge.par, "BGe, order", lambda), results, truegraph)
      
      # GP score, partition
      gp.ofit <- GP.partition.mcmc(data, GP.searchspace, order = F)
      results <- compare_results(gp.ofit, c(gp.par, "GP, partition", lambda), results, truegraph)
      
      # Laplace score, partition
      gp.ofit$weights <- NULL
      gp.ofit$time <- gp.ofit$time - gp.ofit$time2
      results <- compare_results(gp.ofit, c(gp.par, "Laplace, partition", lambda), results, truegraph)
      
      # GP score, order
      gp.pfit <- GP.partition.mcmc(data, GP.searchspace, order = T)
      results <- compare_results(gp.pfit, c(gp.par, "GP, order", lambda), results, truegraph)
      
      # Laplace score, order
      gp.pfit$weights <- NULL
      gp.pfit$time <- gp.pfit$time - gp.pfit$time2
      results <- compare_results(gp.pfit, c(gp.par, "Laplace, order", lambda), results, truegraph)
      
      # DIBS+
      dib.fit <- DiBS(data, par = dib.par)
      results <- compare_results(dib.fit, c(dib.par, "DiBS+", lambda), results, truegraph)
      
      # k-PC, distance correlation
      kPC.dcc.fit <- kPC.dcc.boot(data, nboots = 100, alpha = pc.par)
      results <- compare_results(kPC.dcc.fit, c(pc.par, "kPC-DC", lambda), results, truegraph)
      
      # k-PC, HSIC
      kPC.hsic.fit <- kPC.hsic.boot(data, nboots = 100, alpha = pc.par)
      results <- compare_results(kPC.hsic.fit, c(pc.par, "kPC-HSIC", lambda), results, truegraph)
      
      cat("k =",k,", i =",i,"\n")
    }
  }
  saveRDS(results, paste0("Sims_Results_l=",lambda,".rds"))
}
colnames(results) <- c("ESHD", "eTP", "eFP", "TPR", "FPR_P", 
                       "time", "parameter", "method", "lambda", "graph")
saveRDS(results, "Sims_Results.rds")

# filename <- file.choose()
# results <- readRDS(filename)
