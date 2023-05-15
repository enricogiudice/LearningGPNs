library(tidyverse)
library(ggpubr)

results <- as.data.frame(readRDS(Sims_results.rds))

all_methods <- c("GP, partition", "GP, order", "BGe, partition", "BGe, order", "kPC-HSIC", "kPC-DC", "DiBS+")
some_methods <- c("GP, order", "BGe, order",  "kPC-HSIC", "kPC-DC")
graph_compare <- "dag"  # Type of graph to compare results on 

results %>% filter(method %in% all_methods, graph == graph_compare) %>%
  mutate(ESHD = as.numeric(ESHD),
         TPR = as.numeric(TPR),
         FPRp = as.numeric(FPR_P),
         time = as.numeric(time),
         parameter = as.numeric(parameter),
         lambda = as.numeric(lambda),
         method = factor(method, levels = all_methods)) -> results

lambdas <- c(1, 0.5, 0)  
rocplots <- list()
shdplots <- list()
timeplots <- list()

for(i in 1:length(lambdas)) {
  roc.results <- filter(results, lambda == lambdas[i])
  roc.results %>%
    group_by(method, parameter) %>%
    summarise(TPR = mean(TPR), 
              FPRp = mean(FPRp),
              ESHD = mean(ESHD)) -> summarised_results

  # ROC curve plots
  rocplots[[i]] <- ggplot() + 
    geom_point(data = roc.results %>% 
                filter(method %in% some_methods), 
              aes(x = FPRp, y = TPR, color = method), alpha = 0.15, size = 0.5, shape = 20) + 
    geom_path(data = summarised_results %>% 
                filter(method %in% some_methods), 
              aes(x = FPRp, y = TPR, colour = method), size = 0.6) +
    ggtitle(bquote(lambda ~ "=" ~ .(lambdas[i]))) +
  coord_fixed(ratio = 1, xlim = c(0, 0.45), ylim = c(0.4, 1), expand = FALSE) +  
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank()) +
  scale_color_manual(values=c("#9bb658","#58b6b6","#5866b6","#9b58b6"))

  # ESHD plots
  summarised_results %>%  # Select best tuning parameter for each method
    group_by(method) %>%
    summarise(par = parameter[which.min(ESHD)]) -> par.min
  best.par <- par.min[match(roc.results$method, par.min$method), "par"]
  shd.results <- roc.results[which(roc.results$parameter == best.par), ]
  
  shdplots[[i]] <- ggplot(shd.results, aes(x = method, y = ESHD, color = method)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.3, shape = 16, position = position_jitter(0.1)) +
    scale_color_manual(values=c("#b68158","#9bb658","#58b666","#58b6b6","#5866b6","#9b58b6","#b65881")) +
    ggtitle(bquote(lambda ~ "=" ~ .(lambdas[i]))) +
    xlab("") +
    ylab("E-SHD") +
    theme_light() +
    ylim(c(0, 16)) +  # 30 for n = 15
    theme(legend.position="bottom", legend.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5)) + guides(color = guide_legend(nrow = 1))
  shd_leg <- get_legend(shdplots[[i]])
  shdplots[[i]] <- shdplots[[i]] + theme(legend.position = "none")
  
  # Time plots
  timeplots[[i]] <- ggplot(shd.results, aes(x = method, y = time, color = method)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.3, shape = 16, position = position_jitter(0.1)) +
    scale_color_manual(values=c("#b68158","#9bb658","#58b666","#58b6b6","#5866b6","#9b58b6","#b65881")) +
    ggtitle(bquote(lambda ~ "=" ~ .(lambdas[i]))) +
    xlab("") +
    ylab("Time (s)") +
    theme_light() +
    scale_y_continuous(trans = "log10", limits = c(1, 2200)) +
    theme(legend.position="none", legend.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5)) + guides(color = guide_legend(nrow = 1))
  time_leg <- get_legend(timeplots[[i]])
  timeplots[[i]] <- timeplots[[i]] + theme(legend.position = "none")
}
ggarrange(shdplots[[1]], shdplots[[2]], shdplots[[3]], 
          legend.grob =  shd_leg, ncol = 3, common.legend = T, legend = "bottom")
ggarrange(rocplots[[1]], rocplots[[2]], rocplots[[3]], 
          ncol = 3, common.legend = T, legend = "bottom")
ggarrange(timeplots[[1]], timeplots[[2]], timeplots[[3]], 
          legend.grob = time_leg, ncol = 3, common.legend = T, legend = "bottom")

# size = 9 x 3.6
