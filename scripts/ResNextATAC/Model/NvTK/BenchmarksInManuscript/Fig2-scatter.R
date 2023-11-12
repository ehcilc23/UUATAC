setwd("~/Nvwa_Benchmark/")

library(ggthemes)
library(ggplot2)
library(ggpubr)
library(ggpointdensity) # 绘制密度散点图
library( gridExtra)

######################### scatter smooth ##############################
library(cowplot)
levels_order <- c("mESC_BulkSeq", "mESC_MicrowellSeq", "mESC_DropSeq", 
                  "mESC_MARSseq", "mESC_CELseq", "mESC_SmartSeq","mESC_SCRBseq", "mESC_SmartSeq2", 
                  "mESC_MAGIC-MicrowellSeq",  "mESC_MAGIC-MARSseq", "mESC_MAGIC-DropSeq",
                  "mESC_MAGIC-CELseq", "mESC_MAGIC-SmartSeq", "mESC_MAGIC-SCRBseq", "mESC_MAGIC-SmartSeq2")

generate_p <- function(i){
  print(i)
  prediction <- read.csv(paste0("mESC_benchmark_scTech/", i, "/Test/test_mode_pred_prob.csv"), row.names = 1)
  target_orig <- read.csv(paste0("mESC_benchmark_scTech/", i, "/Test/test_target_prob.csv"), row.names = 1)
  result <- as.data.frame(cbind(target_orig[,1], prediction[,1]))
  colnames(result) <- c("target", "prediction")
  # result <- log(result + 0.001)
  p <- ggplot(result, aes(x=target, y=prediction)) +
    # geom_pointdensity() + #密度散点图（geom_pointdensity）
    geom_point(size=0.5) +
    geom_smooth(method = 'lm') +
    geom_rug() +
    labs(title = i) +#colnames(target_orig)[1]) +
    # theme_base() +
    theme(plot.title = element_text(hjust = 0.5, size = 10))  
  # stat_density2d(colour="grey")
  return(p)
}
ps <- lapply(levels_order, generate_p)


ml <- marrangeGrob(ps, nrow=3, ncol=5)
ml
ggsave("mESC_benchmark_scTech/scatter.all.pdf", ml, width = 12, height = 6.5)