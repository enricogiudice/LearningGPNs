library(kpcalg)

kPC.dcc.boot <- function(data, nboots = 50, alpha = 0.1) {  # k-PC with distance correlation
  mats <- list()  # to store sampled graphs
  N <- nrow(data)
  
  for(i in 1:nboots) {
    boot_data <- data[sample(c(1:N), size = N/2, replace = TRUE), ]
    
    kpc2 <- kpc(suffStat = list(data = boot_data, ic.method = "dcc.perm"),
                indepTest = kernelCItest, alpha = alpha, p = ncol(data))
    mats[[i]] <- as(kpc2@graph, "matrix")
  }
  kpc.main <- kpc(suffStat = list(data = data, ic.method = "dcc.perm"),
                  indepTest = kernelCItest, alpha = alpha, p = ncol(data))
  DAG <- as(kpc.main@graph, "matrix")
  
  return(list(traceadd = list(incidence = mats), DAG = DAG))
}

kPC.hsic.boot <- function(data, nboots = 50, alpha = 0.1) {  # k-PC with hsic
  mats <- list()  
  N <- nrow(data)
  
  for(i in 1:nboots) {
    boot_data <- data[sample(c(1:N), size = N/2, replace = TRUE), ]
    
    kpc2 <- kpc(suffStat = list(data = boot_data, ic.method = "hsic.gamma"),
                indepTest = kernelCItest, alpha = alpha, p = ncol(data))
    mats[[i]] <- as(kpc2@graph, "matrix")
  }
  kpc.main <- kpc(suffStat = list(data = data, ic.method = "hsic.gamma"),
                  indepTest = kernelCItest, alpha = alpha, p = ncol(data))
  DAG <- as(kpc.main@graph, "matrix")
  
  return(list(traceadd = list(incidence = mats), DAG = DAG))
}

DiBS <- function(data, n_particles = 1, p = 1) {  
  n <- ncol(data)
  str_data <- paste(signif(data, 5), collapse = ",")  # convert to strings
  py.args <- paste(str_data, n, n_particles, p)
  
  path <- NULL  # to run on server
  path <- 'cd "/Users/giudic0000/Downloads/Nonlinear scoring/Structure Learning";'
  
  asvec <- system(paste(path, 'python3 algo_dibs.py', py.args), intern = T)
  asnum <- as.numeric(strsplit(asvec, "")[[1]])
  
  if(n_particles == 1) {
    mat = matrix(asnum, ncol = n, byrow = T)
    return(list(traceadd = list(incidence = list(mat)), DAG = mat))
  } 
  else {
    vec_list <- split(asnum, ceiling(seq_along(asnum)/n^2))
    mats <- lapply(vec_list, function(x) matrix(x, ncol = n, byrow = T))
    return(list(traceadd = list(incidence = mats), DAG = mats[[1]]))
  }
}

compareGraphs <- function(fit, truegraph, cpdag = F, absolute = F) {
  DAGList <- fit$traceadd$incidence
  weights <- fit$weights  # normalized log-weights
  SHD <- vector()
  TP <- vector()
  FP <- vector()
  abs <- NULL
  
  for(dag in DAGList) {  # Compute expected SHD, TP & FP
    comp <- compareDAGs(as.matrix(dag), truegraph, cpdag = cpdag)
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
  maxdag <- as.matrix(fit$DAG)
  dagcomp <- compareDAGs(maxdag, truegraph, cpdag)
  
  if(absolute) {
    abs <- c(eTP, eFP)
  }
  
  c(eshd, abs, dagcomp[c("TPR", "FPRn")])
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
