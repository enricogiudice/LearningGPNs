library(gRbase)
library(BiDAG)
library(matrixStats)
library(questionr)
library(ggpubr)
source("Fourier_fns.R")
source("BayesStanFns.R")
source("sampling_fns.R")
source("dualPC.R")
insertSource("GPscore.R", package = "BiDAG")

KL_div <- function(est.post, est.ind, p) {  # reverse KL divergence
  q <- rep(NA, length(p))
  q[est.ind] <- est.post
  d <- q*log(q/p)
  sum(d, na.rm = T)
}

set.seed(101)
lambda <- 1
dual <- T  # use dualPC
n <- 4  # number of nodes
N <- 100  # number of samples

myDAG <- pcalg::randomDAG(n, prob = 0.5, lB = 1, uB = 2) 
trueDAG <- as(myDAG, "matrix")
truegraph <- 1*(trueDAG != 0)
data <- Fou_nldata(truegraph, N, lambda = lambda, noise.sd = 1, standardize = T) 

# List all DAGs with n nodes
all.dags <- list()
adj <- matrix(0, nrow = n, ncol = n)
dag.counter <- 0
all.comb <- rep(list(c(0,1)), n*(n-1))
all.comb <- expand.grid(all.comb)  # all combinations outside of diagonal of adjacency matrix

for(i in 1:nrow(all.comb)) {
  adj[col(adj)!=row(adj)] <- as.numeric(all.comb[i, ])
  
  if(is.DAG(adj)) {
    dag.counter <- dag.counter + 1
    all.dags[[dag.counter]] <- adj
  }
}

# Compute true posterior (and save all scores)
true.post <- rep(NA, dag.counter)
parent.scores <- data.frame(sets = character(0), newscore = character(0))  
parent.scores <- rep(list(parent.scores), n)

for(k in 1:dag.counter) {
  dag <- all.dags[[k]]
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
  true.post[k] <- curr_score
}

# Order true posterior
true.order <- order(true.post, decreasing = T)
true.post <- true.post[true.order]
true.p <- exp(true.post - logSumExp(true.post))
all.dags <- all.dags[true.order]

# Estimate posterior via order GP sampling
toburn <- 250
GP.searchspace <- set.searchspace(data, dual = T, "GP", alpha = 0.3)
iters <- c(30e1, 44e1, 65e1, 96e1, 14e2, 21e2, 31e2, 46e2, 67e2, 10e3)  # could do more...
results <- data.frame()

for(i in 1:length(iters)) {
  parfit <- partitionMCMC(GP.searchspace$score, alpha = 0.2, startDAG = GP.searchspace$DAG,
                      scoretable = GP.searchspace$scoretable, startspace = GP.searchspace$endspace, 
                      iterations = 2 * (iters[i] + toburn), stepsave = 2)
  sampled_dags <- parfit$traceadd$incidence[-(1:toburn)]

  # Find index in all.dags of sampled dags 
  all.vecdags <- lapply(all.dags, c)
  post.indexes <- sapply(sampled_dags, function(x) which(all.vecdags %in% list(as.numeric(x))))

  # Compute normalized weights for sampled DAGs
  ord_trace <- parfit$trace[-(1:toburn)]
  weights <- true.post[post.indexes] - ord_trace
  weights <- weights - logSumExp(weights)

  # Normalize posteriors
  tab.ind <- wtd.table(post.indexes, weights = exp(weights))
  est.post <- as.numeric(tab.ind)
  est.ind <- as.numeric(names(tab.ind))

  # Laplace approximate posterior
  lap.tab <- table(post.indexes)
  lap.post <- as.numeric(lap.tab)/sum(lap.tab)
  
  # Save results
  results <- rbind(results, data.frame(kl = KL_div(est.post, est.ind, true.p), 
                                       iter = iters[i], group = "GP, weighted"))
  results <- rbind(results, data.frame(kl = KL_div(lap.post, est.ind, true.p), 
                                       iter = iters[i], group = "GP, Laplace"))
}

# Plot KL divergences and posterior for last setting
post_plots <- list()
post_plots[[1]] <- ggplot(results, aes(x = iter, y = kl, group = group)) +
  geom_line(aes(linetype = group)) +
  geom_point(aes(shape = group, color = group)) +
  scale_shape_manual(values = c(2, 3)) +
  scale_color_manual(values = c('#00c700','#db0000')) +
  scale_linetype_manual(values = c("twodash", "dashed")) +
  scale_x_continuous(trans = 'log10') +
  xlab("Samples") + ylab("K-L divergence") +
  theme_light() +
  theme(legend.title = element_blank(), legend.position = "none")

post_data <- rbind(data.frame(x = est.ind, y = est.post, group = "GP, weighted"),
                   data.frame(x = est.ind, y = lap.post, group = "GP, Laplace"))
post_plots[[2]] <- ggplot() +
  geom_bar(data = data.frame(x = 1:dag.counter, y = true.p), aes(x,y), 
           stat = 'identity', width = 0.6) +
  geom_point(data = post_data, aes(x, y, shape = group, color = group), size = 1) +
  scale_shape_manual(values = c(2, 3)) +
  scale_color_manual(values = c('#00c700','#db0000')) +
  xlab("DAG") + ylab("Posterior probability") +
  theme_light() +
  theme(legend.title = element_blank(), legend.position = "none") +
  scale_y_sqrt()

ggarrange(post_plots[[1]], post_plots[[2]], ncol = 2, widths = c(1, 2), 
          common.legend = T, legend = "bottom")
# size = 3.3 x 9
