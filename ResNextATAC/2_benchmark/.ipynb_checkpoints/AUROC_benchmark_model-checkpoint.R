###### set working dir ######
parser <- ArgumentParser()
parser$add_argument("-w", "--workdir", required=TRUE, help="absolute path of working directory")
args <- parser$parse_args()
workdir <- args$workdir
setwd(workdir)


library(RColorBrewer)
library(ggplot2)

###### compare auroc results between models ######
Metric1 <- read.csv('./resnext_Metric.csv',row.names = 1)
Metric1 <- as.data.frame(t(as.matrix(Metric1)))
Metric1$group <- 'ResNext'

Metric2 <- read.csv('./scbasset_Metric.csv',row.names = 1)
Metric2 <- as.data.frame(t(as.matrix(Metric2)))
Metric2$group <- 'scBasset'

Metric3 <- read.csv('./nvwa_Metric.csv',row.names = 1)
Metric3 <- as.data.frame(t(as.matrix(Metric3)))
Metric3$group <- 'nvwa'


###### AUROC plot ######
Metric_total <- rbind(Metric1,Metric2,Metric3)

Metric_total$group=factor(Metric_total$group,levels = c('ResNext','nvwa','scBasset'))
head(Metric_total)
p1 <- ggplot(Metric_total, aes(x=group, y=auroc,color=group)) + 
  geom_boxplot(outlier.size = 0.5)+
  scale_color_manual(values= brewer.pal(11,'Spectral')[8:10]) +
  theme(axis.text = element_text(size= 11, colour = "black"),
        axis.title = element_text(size = 14, colour = "black", face = "bold"), # coordinate aesthetic
        plot.title = element_text(size = 20, hjust = 0.5), # main title aesthetic
        # panel.background = element_blank(),
        # panel.grid = element_blank(),
        plot.margin = unit(rep(2, 4), 'lines')) +
  labs(y = "AUROC", x = "",title = '')+geom_hline(yintercept =0.8 ,lty='dashed',sie=0.25)+ylim(0.5,1)
p1
pdf(file = paste('AUROC_benchmark_model.pdf'),width = 10,height = 8)
print(p1)
dev.off()
