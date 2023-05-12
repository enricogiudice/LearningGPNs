# This script looks at the behavior of the GP score across different levels of non-Gaussianity
library(tidyverse)
library(BiDAG)

source("Fourier_fns.R")
source("BayesStanFns.R")
source("sampling_fns.R")
source("dualPC.R")
insertSource("GPscore.R", package = "BiDAG")

# Compute score of DAG
GPscoreDAG <- function(dag, data) {
  n <- ncol(data)
  curr_score <- 0
  
  for(x in 1:n) {
    pax <- dag[,x]
    
    if(sum(pax) == 0) {  # No parents
      loc_score <- Gauss.mcmc(data[,x])
    }
    
    else {  # One or more parents
      loc_score <- GP.mcmc(data[ ,x], data[ ,which(pax==1)])
    }
    curr_score <- curr_score + loc_score  # build score
  }
  return(curr_score)
}

init.seed <- 100
lambdas <- seq(0, 1, length.out = 19)  # non-linearity; zero is linear
iter <- 100  # number of simulations per value of lambda
dual <- T  # use dualPC
n <- 5  # number of nodes
N <- 100  # number of samples
results <- data.frame()

for(l in lambdas) {
  
  for (i in 1:iter) {
    set.seed(init.seed+i)

    # Generate DAG & data
    truegraph <- rbind(cbind(rep(0, n-1), diag(n-1)), rep(0, n))
    reversegraph <- t(truegraph)
    data <- Fou_nldata(truegraph, N, lambda = l, noise.sd = 1, standardize = T) 
  
    # GP, partition
    score <- scoreparameters("usr", data, usrpar = list(pctesttype = "bge"))
    scorediff <- GPscoreDAG(truegraph, data) - GPscoreDAG(reversegraph, data)
  
    results <- rbind(results, c(scorediff, l))
  }
}
saveRDS(results, "Scorequi_results.rds")


# Plot results
colnames(results) <- c("Scorediff", "lambda")
results %>%
  group_by(lambda) %>%
  summarise(mediff = median(Scorediff),
            hi = quantile(Scorediff, 0.9),
            lo = quantile(Scorediff, 0.1)) -> avgresults

ggplot(avgresults, aes(x = lambda, y = mediff)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size=0.7, color="#af3c3c", fill="white", shape=16) +
  xlab(bquote(lambda)) + ylab("Difference in scores") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid.minor = element_blank(),
        panel.border = element_blank())

# size: 3 x 4.1
