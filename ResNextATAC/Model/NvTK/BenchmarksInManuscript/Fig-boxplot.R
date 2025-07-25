setwd("~/Nvwa_Benchmark/Benchmark_scTech_mESC/")

library(ggthemes)
library(ggplot2)
library(ggpubr)

levels_order <- c("mESC_BulkSeq", "mESC_MicrowellSeq", "mESC_DropSeq", 
            "mESC_MARSseq", "mESC_CELseq", "mESC_SmartSeq", "mESC_SCRBseq", "mESC_SmartSeq2")

# levels_order <- c("mESC_CELseq2A", "mESC_CELseq2B", "mESC_DropSeqA", "mESC_DropSeqB", "mESC_MARSseqA", "mESC_MARSseqB",
#                   "mESC_SCRBseqA", "mESC_SCRBseqB", "mESC_SmartSeq2A", "mESC_SmartSeq2B", "mESC_SmartSeqA", "mESC_SmartSeqB") 
#                   
# levels_order <- c("mESC_BulkSeq", "mESC_MAGIC-MicrowellSeq",  "mESC_MAGIC-MARSseq", "mESC_MAGIC-DropSeq",
#             "mESC_MAGIC-CELseq", "mESC_MAGIC-SmartSeq", "mESC_MAGIC-SCRBseq", "mESC_MAGIC-SmartSeq2")


###########################################
metric.fpaths <- paste0(levels_order, "/Test/Metric.csv")
metric.fpaths <- metric.fpaths[file.exists(metric.fpaths)]
metric.list <- lapply(metric.fpaths, function(metric.fname){
  print(metric.fname)
  metric <- as.data.frame(t(read.csv(metric.fname, row.names = 1)))
  metric$type <- strsplit(metric.fname, '/')[[1]][1]
  return(metric)
})
metric <- do.call(rbind, metric.list)
# metric$type <- gsub("_Trial.", "", metric$type)
unique(metric$type)

# sort factor
metric$type <- factor(metric$type, levels = levels_order)
metric <- metric[order(metric$type),]
# metric$value <- ifelse(metric$value > 0.5, metric$value, 1 - metric$value)
dim(metric); head(metric)

pdf("./PCC_mESC_benchmark_scTech.pdf", width = 10, height = 6)
ggplot(metric, aes(x=type, y=pcc, fill=type)) + 
  geom_boxplot() + theme_base() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=10), legend.position="none") +
  ylab("PCC") #+
# stat_compare_means(aes(group = type), method="wilcox.test", label = "p.format", label.y = 1) 
dev.off()

pdf("./SCC_mESC_benchmark_scTech.pdf", width = 10, height = 6)
ggplot(metric, aes(x=type, y=pcc, fill=type)) + 
  geom_boxplot() + theme_base() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=10), legend.position="none") +
  ylab("SCC") #+
# stat_compare_means(aes(group = type), method="wilcox.test", label = "p.format", label.y = 1) 
dev.off()

###########################################
IC.fpaths <- paste0(levels_order, "/Motif/meme_IC.csv")
IC.fpaths <- IC.fpaths[file.exists(IC.fpaths)]
IC.list <- lapply(IC.fpaths, function(IC.fname){
  print(IC.fname)
  IC <- read.csv(IC.fname)
  IC$type <- strsplit(IC.fname, '/')[[1]][1]
  return(IC)
})
IC <- do.call(rbind, IC.list)
# IC$type <- gsub("_Trial.", "", IC$type)
unique(IC$type)

# sort factor
IC$type <- factor(IC$type, levels =  levels_order)
IC <- IC[order(IC$type),]
# IC$IC <- ifelse(IC$IC > 30, 30, IC$IC)
# IC$IC <- ifelse(IC$IC < 5, 5, IC$IC)
dim(IC); head(IC)

pdf("./IC_mESC_benchmark_scTech.pdf", width = 10, height = 6)
ggplot(IC, aes(x=type, y=IC, fill=type)) + 
  geom_boxplot() + theme_base() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=10), legend.position="none") +
  ylab("IC(Information Content)") #+
# stat_compare_means(aes(group = type), method="wilcox.test", label = "p.format", label.y = 1) 
dev.off()