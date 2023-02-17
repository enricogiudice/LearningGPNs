usrDAGcorescore <- function (j, parentnodes, n, param) {  # GP Laplace approximate score 
  lp <- length(parentnodes)  # number of parents
  scoreconstvec <- param$scoreconstvec
  y <- param$data[ ,j]
  switch(as.character(lp),
         "0" = {  # no parents
           llik <- Gauss.lap(y)
           corescore <- scoreconstvec[lp+1] + llik
         },
         {  # one or more parents
           corescore <- scoreconstvec[lp+1]
           X <- param$data[ ,parentnodes]
           llik <- GP.lap(y, X)
           corescore <- corescore + llik
           if (!is.null(param$logedgepmat)) {  # if there is an additional edge penalization
             corescore <- corescore - sum(param$logedgepmat[parentnodes, j])
           }
         })
  return(corescore)
}
