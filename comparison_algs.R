library(kpcalg)

kPC.dcc.boot <- function(data, nboots = 100, alpha = 0.1) {  # k-PC with distance correlation
  start <- Sys.time()
  mats <- list()  # to store sampled graphs
  N <- nrow(data)
  
  for(i in 1:nboots) {
    boot_data <- data[sample(c(1:N), size = N/2, replace = TRUE), ]
    
    kpc2 <- kpc(suffStat = list(data = boot_data, ic.method = "dcc.perm"), 
                m.max = 4, indepTest = kernelCItest, alpha = alpha, p = ncol(data))
    mats[[i]] <- as(kpc2@graph, "matrix")
  }
  kpc.main <- kpc(suffStat = list(data = data, ic.method = "dcc.perm"),
                  indepTest = kernelCItest, alpha = alpha, p = ncol(data))
  DAG <- as(kpc.main@graph, "matrix")
  time <- Sys.time() - start
  
  return(list(traceadd = list(incidence = mats), DAG = DAG, 
              time = as.numeric(time, units = "secs")))
}

kPC.hsic.boot <- function(data, nboots = 100, alpha = 0.1) {  # k-PC with hsic
  start <- Sys.time()
  mats <- list()  
  N <- nrow(data)
  
  for(i in 1:nboots) {
    boot_data <- data[sample(c(1:N), size = N/2, replace = TRUE), ]
    
    kpc2 <- kpc(suffStat = list(data = boot_data, ic.method = "hsic.gamma"), 
                m.max = 4, indepTest = kernelCItest, alpha = alpha, p = ncol(data))
    mats[[i]] <- as(kpc2@graph, "matrix")
  }
  kpc.main <- kpc(suffStat = list(data = data, ic.method = "hsic.gamma"),
                  indepTest = kernelCItest, alpha = alpha, p = ncol(data))
  DAG <- as(kpc.main@graph, "matrix")
  time <- Sys.time() - start
  
  return(list(traceadd = list(incidence = mats), DAG = DAG, 
              time = as.numeric(time, units = "secs")))
}

DiBS <- function(data, n_particles = 1, par = 1) {
  start <- Sys.time()
  n <- ncol(data)
  str_data <- paste(signif(data, 5), collapse = ",")  # convert to strings
  py.args <- paste(str_data, n, n_particles, par)
  
  path <- NULL  # to run on server
  path <- 'cd "/Users/giudic0000/Downloads/Nonlinear scoring/Structure Learning";'
  
  asvec <- system(paste(path, 'python3 algo_dibs.py', py.args), intern = T)
  asnum <- as.numeric(strsplit(asvec, "")[[1]])
  time <- Sys.time() - start
  
  if(n_particles == 1) {
    mat = matrix(asnum, ncol = n, byrow = T)
    return(list(traceadd = list(incidence = list(mat)), DAG = mat, 
                time = as.numeric(time, units = "secs")))
  } 
  else {
    vec_list <- split(asnum, ceiling(seq_along(asnum)/n^2))
    mats <- lapply(vec_list, function(x) matrix(x, ncol = n, byrow = T))
    return(list(traceadd = list(incidence = mats), DAG = mats[[1]], 
                time = as.numeric(time, units = "secs")))
  }
}

# Compare estimated cpdag, dag, skeleton & pattern graph to truth
compare_results <- function(fit, method_vec, result_df, trueDAG) {
  new <- length(method_vec)
  
  method_vec[new+1] <- "dag"  # Compare dag
  result_dag <- compareFit(fit$traceadd$incidence, fit$DAG, fit$time, trueDAG, fit$weights)
  result_df <- rbind(result_df, data.frame(t(c(result_dag, method_vec))))
  
  method_vec[new+1] <- "cpdag"  # Compare cpdag
  cpdagList <- NULL  # check if DAGs are acyclic
  try(cpdagList <- lapply(fit$traceadd$incidence, function(x) BiDAG:::dagadj2cpadj(x)), silent = T)
  if(!is.null(cpdagList)) {    # check if DAGs are acyclic
    trueCPDAG <- BiDAG:::dagadj2cpadj(trueDAG)
    result_cpdag <- compareFit(cpdagList, BiDAG:::dagadj2cpadj(fit$DAG),
                              fit$time, trueCPDAG, fit$weights)
    result_df <- rbind(result_df, data.frame(t(c(result_cpdag, method_vec))))
  
    method_vec[new+1] <- "skeleton"  # Compare skeleton
    skelList <- lapply(fit$traceadd$incidence, function(x) Gskel(as.matrix(x)))
    trueskel <- Gskel(trueDAG)
    result_skel <- compareFit(skelList, Gskel(as.matrix(fit$DAG)), fit$time, trueskel, fit$weights)
    result_df <- rbind(result_df, data.frame(t(c(result_skel, method_vec))))
    
    method_vec[new+1] <- "pattern"  # Compare pattern graph
    pattList <- lapply(cpdagList, function(x) pdag2pattern(x))
    truepatt <- pdag2pattern(trueCPDAG)
    result_patt <- compareFit(pattList, pdag2pattern(as.matrix(fit$DAG)), fit$time, truepatt, fit$weights)
    result_df <- rbind(result_df, data.frame(t(c(result_patt, method_vec))))
  }
  return(result_df)
}

# Compare the fit to the truth
compareFit <- function(graphlist, MAPgraph, time, truegraph, weights = NULL) {
  SHD <- vector()
  TP <- vector()
  FP <- vector()
  
  for(graph in graphlist) {  # Compute expected SHD, TP & FP
    comp <- compareGs(as.matrix(graph), truegraph)
    SHD <- append(SHD, comp["SHD"])
    TP <- append(TP, comp["TP"])
    FP <- append(FP, comp["FP"])
  }
  
  if(!is.null(weights)) {
    eshd <- sum(SHD * exp(weights))
    eTP <- sum(TP * exp(weights))
    eFP <- sum(FP * exp(weights))
  }
  
  else {
    eshd <- mean(SHD)
    eTP <- mean(TP)
    eFP <- mean(FP)
  }
  maxdag <- as.matrix(MAPgraph)
  graphcomp <- compareGs(maxdag, truegraph)
  
  c(eshd, eTP, eFP, graphcomp[c("TPR", "FPR_P")], time)
}

post.edges <- function(fit, weights = NULL) {
  MCMCchain <- fit$traceadd$incidence
  endstep <- length(MCMCchain)
  
  if(is.null(weights)) {
    w <- rep(1/endstep, endstep)
  }
  
  else {
    w <- exp(weights)
  }
  wchain <- Map('*', MCMCchain, as.list(w))
  incidence <- as.matrix(Reduce("+", wchain[1:endstep]))
  
  c(incidence[6,7], 1-incidence[6,8])  # probability of edges Erk -> Akt and Erk -|> PKA
}

# This function extracts the skeleton from a graph
Gskel <- function(incidence) {
  1*(incidence|t(incidence))
}

# This function compares an estimated graph to the true one
compareGs <- function (estG, trueG) {
  estSkel <- Gskel(estG) # estimated skeleton
  trueSkel <- Gskel(trueG) # true skeleton
  P <- sum(trueSkel)/2 # number of positives
  diffSkel <- estSkel - trueSkel
  extra_edges <- which(diffSkel > 0) # edges in estimated but not true EG
  FP <- length(extra_edges)/2 # count to FPs
  estG[extra_edges] <- 0 # remove them from further comparisons
  missing_edges <- which(diffSkel < 0) # edges in true but not estimated EG
  FN <- length(missing_edges)/2 # count to FNs
  trueG[missing_edges] <- 0 # remove them from further comparisons
  # modified graphs have the same skeletons, so now just need to count mismatches
  mismatches <- 1*(estG != trueG)
  wrong_order <- sum(Gskel(mismatches))/2 # number of wrongly oriented edges
  FP <- FP + wrong_order/2 # include half in FP
  FN <- FN + wrong_order/2 # and half in FN
  SHD <- FP + FN # shd is the sum of errors
  TP <- P - FN # true positives are without false negatives
  # TPR, FPR_P
  if (P == 0) { # true graph is empty
    if (FP >= 0) {
      TPR <- 0
      FPR_P <- 1
    } else {
      TPR <- 1
      FPR_P <- 0
    }
  } else { # true graph is non-empty
    TPR <- TP/P
    FPR_P <- FP/P
  }
  compGs <- c(TP, FP, SHD, TPR, FPR_P, P)
  names(compGs) <- c("TP","FP", "SHD", "TPR", "FPR_P", "P")
  return(compGs)
}
