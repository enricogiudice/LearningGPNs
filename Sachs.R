library(BiDAG)
library(matrixStats)

source("/Users/giudic0000/Downloads/Nonlinear scoring/Structure Learning/BayesStanFns.R")
source("/Users/giudic0000/Downloads/Nonlinear scoring/Structure Learning/sampling_fns.R")
source("/Users/giudic0000/Downloads/Nonlinear scoring/Structure Learning/comparison_algs.R")
source("/Users/giudic0000/Downloads/Nonlinear scoring/Structure Learning/dualPC.R")
insertSource("~/Downloads/Nonlinear scoring/Structure Learning/GPscore.R", package = "BiDAG")
sachs.all <- read.delim("~/Downloads/Nonlinear scoring/Structure Learning/sachs_allexps.csv")
trueDAGbn <- readRDS("~/Downloads/Nonlinear scoring/Structure Learning/sachs_graph.rds")
set.seed(100)
results <- data.frame()
exp_num <- 1

sachs.exp <- subset(sachs.all[sachs.all$experiment == exp_num, ], select = -experiment) 
sachs.data <- scale(log2(sachs.exp + 0.5))

# Set initial search spaces
GP.searchspace = set.searchspace(sachs.data, dual = T, "GP")  # more correct
bge.searchspace = set.searchspace(sachs.data, dual = T, "bge")

# BGe score, partition
bge.fit <-  bge.partition.mcmc(bge.searchspace, order = F, iterations = 1200)
bge.edgep <- post.edges(bge.fit)
results <- compare_results(bge.fit, c(bge.edgep, "BGe, partition"), results, trueDAGbn)

# BGe score, order
bge.fit <-  bge.partition.mcmc(bge.searchspace, order = T, iterations = 1200)
bge.edgep <- post.edges(bge.fit)
results <- compare_results(bge.fit, c(bge.edgep, "BGe, order"), results, trueDAGbn)

# GP score, partiton
gp.fit <-  GP.partition.mcmc(sachs.data, GP.searchspace, order = F, iterations = 1200)
gp.edgep <- post.edges(gp.fit, gp.fit$weights)
results <- compare_results(gp.fit, c(gp.edgep, "GP, partition"), results, trueDAGbn)

# GP score, order
gp.fit <-  GP.partition.mcmc(sachs.data, GP.searchspace, order = T, iterations = 1200)
gp.edgep <- post.edges(gp.fit, gp.fit$weights)
results <- compare_results(gp.fit, c(gp.edgep, "GP, order"), results, trueDAGbn)

# GP with interactions, partition
source("/Users/giudic0000/Downloads/Nonlinear scoring/Structure Learning/BayesStanFns_interact.R")
gp_int.fit <- GP.partition.inter(sachs.data, dual = T, order = F)
gp_int.edgep <- post.edges(gp_int.fit, gp_int.fit$weights)
results <- compare_results(gp_int.fit, c(gp_int.edgep, "GP_int, partition"), results, trueDAGbn)

# GP with interactions, order
gp_int.fit <- GP.partition.inter(sachs.data, dual = T, order = T)
gp_int.edgep <- post.edges(gp_int.fit, gp_int.fit$weights)
results <- compare_results(gp_int.fit, c(gp_int.edgep, "GP_int, order"), results, trueDAGbn)

# DIBS+
dib.fit <- DiBS(as.matrix(sachs.data), 20)
dib.edgep <- post.edges(dib.fit)
results <- compare_results(dib.fit, c(dib.edgep, "DiBS+"), results, trueDAGbn)

# k-PC, distance correlation
kPC.dcc.fit <- kPC.dcc.boot(sachs.data, nboots = 200)
kPC.dcc.edgep <- post.edges(kPC.dcc.fit)
results <- compare_results(kPC.dcc.fit, c(kPC.dcc.edgep, "kPC-DC"), results, trueDAGbn)

# k-PC, HSIC
kPC.hsic.fit <- kPC.hsic.boot(sachs.data, nboots = 200)
kPC.hsic.edgep <- post.edges(kPC.hsic.fit)
results <- compare_results(kPC.hsic.fit, c(kPC.hsic.edgep, "kPC-HSIC"), results, trueDAGbn)

colnames(results) <- c("ESHD", "eTP", "eFP", "TPR", "FPR_P", "time",
                       "ErktoAkt", "ErktoPKA", "Scorefn", "graph")
saveRDS(results, "Sachs_results.rds")
