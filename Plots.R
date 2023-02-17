library(tidyverse)
library(ggpubr)

readRDS("Results/Sims_Results.rds")

all_methods <- c("GP, partition", "GP, order", "BGe, partition", "BGe, order", "DiBS+", "kPC-HSIC", "kPC-DC")
some_methods <- c("GP, order", "BGe, order",  "kPC-HSIC", "kPC-DC")
results  <- mutate(results, 
                   ESHD = as.numeric(ESHD),
                   Lambda = as.numeric(Lambda),
                   FPRp = as.numeric(results$FPRn),
                   TPR = as.numeric(results$TPR),
                   parameter = as.numeric(parameter),
                   method = factor(Scorefn, levels = all_methods))

lambdas <- c(0, 0.5, 1)  
rocplots <- list()
shdplots <- list()

for(i in 1:length(lambdas)) {
  roc.results <- filter(results, Lambda == lambdas[i])
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
              aes(x = FPRp, y = TPR, colour = method), size = 0.5) +
    ggtitle(bquote(lambda ~ "=" ~ .(lambdas[i]))) +
  coord_fixed(ratio = 1, xlim = c(0, 0.45), ylim = c(0.45, 1), expand = FALSE) +  # try unfixing the ratio
    theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("#b69f58", "#58b678", "#5878b6", "#b358b6"))
  
  
  # ESHD plots
  summarised_results %>%  # Select best tuning parameter for each method
    group_by(method) %>%
    summarise(par = parameter[which.min(ESHD)]) -> par.min
  best.par <- par.min[match(roc.results$method, par.min$method), "par"]
  shd.results <- roc.results[which(roc.results$parameter == best.par), ]
  
  shdplots[[i]] <- ggplot(shd.results, aes(x = method, y = ESHD, color = method)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.3, shape = 16, position = position_jitter(0.1)) +
    scale_color_manual(values=c("#b65858", "#b69f58", "#89b658", "#58b678", "#58adb6", "#5878b6", "#b358b6")) +
    ggtitle(bquote(lambda ~ "=" ~ .(lambdas[i]))) +
    xlab("") +
    ylab("E-SHD") +
    theme_light() +
    ylim(c(0, 16)) +
    theme(legend.position="none", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(hjust = 0.5)) 
}
ggarrange(shdplots[[1]], shdplots[[2]], shdplots[[3]], 
          ncol = 3, common.legend = T, legend = "bottom")
ggarrange(rocplots[[1]], rocplots[[2]], rocplots[[3]], ncol = 3, common.legend = T, legend = "bottom")
