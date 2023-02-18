library(BiDAG)
library(matrixStats)

source("BayesStanFns.R")
source("sampling_fns.R")
source("comparison_algs.R")
source("dualPC.R")
insertSource("GPscore.R", package = "BiDAG")
sachs.all <- read.delim("sachs_allexps.csv")
trueDAGbn <- readRDS("sachs_graph.rds")
set.seed(101)
results <- data.frame()
exp_num <- 1

sachs.exp <- subset(sachs.all[sachs.all$experiment == exp_num, ], select = -experiment) 
sachs.data <- scale(log2(sachs.exp + 0.5))

# Set initial search spaces
GP.searchspace = list(score = scoreparameters("usr", sachs.data, 
                                              usrpar = list(pctesttype = "bge")))
bge.searchspace = list(score = scoreparameters("bge", sachs.data))

# BGe score, partition
bge.fit <-  bge.partition.mcmc(bge.searchspace, order = F, iterations = 2000)
bge.comp <- compareGraphs(bge.fit, trueDAGbn, absolute = T)
bge.edgep <- post.edges(bge.fit)
results <- rbind(results, c(as.character(bge.comp), bge.edgep, "BGe, partition", exp_num))

# BGe score, order
bge.fit <-  bge.partition.mcmc(bge.searchspace, order = T, iterations = 2000)
bge.comp <- compareGraphs(bge.fit, trueDAGbn, absolute = T)
bge.edgep <- post.edges(bge.fit)
results <- rbind(results, c(as.character(bge.comp), bge.edgep, "BGe, order", exp_num))
  
# GP score, partiton
gp.fit <-  GP.partition.mcmc(sachs.data, GP.searchspace, order = F, iterations = 2000)
gp.comp <- compareGraphs(gp.fit, trueDAGbn, absolute = T)
gp.edgep <- post.edges(gp.fit, gp.fit$weights)
results <- rbind(results, c(as.character(gp.comp), gp.edgep, "GP, partition", exp_num))

# GP score, order
gp.fit <-  GP.partition.mcmc(sachs.data, GP.searchspace, order = T, iterations = 2000)
gp.comp <- compareGraphs(gp.fit, trueDAGbn, absolute = T)
gp.edgep <- post.edges(gp.fit, gp.fit$weights)
results <- rbind(results, c(as.character(gp.comp), gp.edgep, "GP, order", exp_num))

# DIBS+
dib.fit.DAG <- DiBS(as.matrix(sachs.data), 1)
dib.fit <- DiBS(as.matrix(sachs.data), 20)
dib.comp <- compareGraphs(dib.fit, trueDAGbn, absolute = T)
dib.edgep <- post.edges(dib.fit)
results <- rbind(results, c(as.character(dib.comp), dib.edgep, "DiBS+", exp_num))

# k-PC, distance correlation
kPC.dcc.fit <- kPC.dcc.boot(sachs.data, nboots = 200)
kPC.dcc.comp <- compareGraphs(kPC.dcc.fit, trueDAGbn, absolute = T)
kPC.dcc.edgep <- post.edges(kPC.dcc.fit)
results <- rbind(results, c(as.character(kPC.dcc.comp), kPC.dcc.edgep, "kPC-DC", exp_num))

# k-PC, HSIC
kPC.hsic.fit <- kPC.hsic.boot(sachs.data, nboots = 200)
kPC.hsic.comp <- compareGraphs(kPC.hsic.fit, trueDAGbn, absolute = T)
kPC.hsic.edgep <- post.edges(kPC.hsic.fit)
results <- rbind(results, c(as.character(kPC.hsic.comp), kPC.hsic.edgep, "kPC-HSIC", exp_num))

colnames(results) <- c("ESHD", "ETP", "EFP", "TPR(MAP)", "FPRp(MAP)",
                       "ErktoAkt", "ErktoPKA", "Scorefn", "Experiment")
saveRDS(results, "Results/Sachs_results.rds")
