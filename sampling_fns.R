set.searchspace <- function(data, dual, method, par = 1, alpha = 0.05) {
  start <- Sys.time()
  startspace <- NULL
  
  if(dual) {
    cor_mat <- cor(data)
    startspace <- dual_pc(cor_mat, nrow(data), alpha = alpha, skeleton = T)
  }
  
  if(method == "GP") {
    n <- ncol(data)
    score <- scoreparameters("usr", data, 
                             usrpar = list(pctesttype = "bge", edgepmat = matrix(par, n, n)))
  }
  
  if(method == "bge") {
    score <- scoreparameters("bge", data, bgepar = list(am = par))
  }
  searchspace <- iterativeMCMC(scorepar = score, startspace = startspace, hardlimit = 14, 
                               verbose = F, scoreout = TRUE, alphainit = 0.01)
  time <- Sys.time() - start
  
  list(score = score, scoretable = searchspace$scoretable, DAG = searchspace$DAG, 
       maxorder = searchspace$maxorder, endspace = searchspace$endspace, time = time)
}

GP.partition.mcmc <- function(data, searchspace, alpha = 0.05, 
                              order = FALSE, burnin = 0.33, iterations = 600) {
  start <- Sys.time()
  n <- ncol(data)
  myScore <- searchspace$score
  if(order) {
    fit <- orderMCMC(myScore, MAP = FALSE, chainout = TRUE, alpha = alpha, 
                     startorder = searchspace$maxorder, scoretable = searchspace$scoretable, 
                     startspace = searchspace$endspace, iterations = iterations, stepsave = 4)
  }
  else {
    fit <- partitionMCMC(myScore, alpha = alpha, startDAG = searchspace$DAG,
                        scoretable = searchspace$scoretable, startspace = searchspace$endspace,
                        iterations = iterations, stepsave = 4)
  }
  inter <- Sys.time()
  # Make list of sampled parents for each variable
  parent.scores <- data.frame(sets = character(0), newscore = character(0))  
  parent.scores <- rep(list(parent.scores), n)
  weights <- NULL  # store importance weights for all sampled DAGs
  toburn <- round(burnin * fit$info$samplesteps)
  fit$traceadd$incidence <- fit$traceadd$incidence[-(1:toburn)]
  fit$trace <- fit$trace[-(1:toburn)]
  ndags <- length(fit$trace)
  
  for(k in 1:ndags) {
    dag <- fit$traceadd$incidence[[k]]
    curr_score <- 0
    
    for(x in 1:n) {
      set <- parent.scores[[x]]$sets
      pax <- dag[,x]
      pax.str <- paste(pax, collapse = "")  # parents of x
      check <- is.na(match(set, pax.str))  # check if parent set is already there
      
      if(all(check)) {  # new parent set
        
        if(sum(pax) == 0) {  # Compute local score
            loc_score <- Gauss.mcmc(data[,x])
        }
        else {
          loc_score <- GP.mcmc(data[ ,x], data[ ,which(pax==1)])
        }
        parent.scores[[x]][length(set)+1, ] <- c(pax.str, loc_score)
      }
      
      else {  # fetch score from parent.scores
        ind <- which(check == F)  # index of pax in set
        loc_score <- as.numeric(parent.scores[[x]]$newscore[ind])
      }
      curr_score <- curr_score + loc_score  # build score
    }
    weights[k] <- curr_score - fit$trace[k]   # weight for current DAG 
  }
  
  fit$weights <- weights - logSumExp(weights)  # normalize weights
  end <- Sys.time()
  time2 <- end - inter
  time <- end - start + searchspace$time
  fit$time <- as.numeric(time, units = "secs")
  fit$time2 <- as.numeric(time2, units = "secs")
  
  return(fit)
}

bge.partition.mcmc <- function(searchspace, alpha = 0.05, 
                               order = FALSE, burnin = 0.33, iterations = 600) {
  start <- Sys.time()
  BGEScore <- searchspace$score
  
  if(order) {
    bge.fit <- orderMCMC(BGEScore, MAP = FALSE, chainout = TRUE, alpha = alpha, 
                         startorder = searchspace$maxorder, scoretable = searchspace$scoretable,
                         startspace = searchspace$endspace, iterations = iterations, stepsave = 4)
  }
  else {
    bge.fit <- partitionMCMC(BGEScore, alpha = alpha, startDAG = searchspace$DAG, 
                             scoretable = searchspace$scoretable, startspace = searchspace$endspace,
                             iterations = iterations, stepsave = 4)
  }
  toburn <- round(burnin * bge.fit$info$samplesteps)
  
  bge.fit$traceadd$incidence <- bge.fit$traceadd$incidence[-(1:toburn)]
  time <- Sys.time() - start + searchspace$time
  bge.fit$time <- as.numeric(time, units = "secs")
  
  return(bge.fit)
}

GP.partition.inter <- function(data, dual, order = F, par = 1) {
  GP.searchspace = set.searchspace(data, dual, "GP", par = par)
  gp.fit <-  GP.partition.mcmc(data, GP.searchspace, order = order, iterations = 1200)

  return(gp.fit)
}
