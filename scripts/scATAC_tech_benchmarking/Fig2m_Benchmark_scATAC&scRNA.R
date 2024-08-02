library(data.table)
library(dplyr)

### scATAC
df_UU <- fread('metadata_UUATAC.txt') %>% as.data.frame()
stat_UU <- data.frame(nMolecule = sum(df_UU$UQ), dataset = 'UUATAC')

df_CH <- fread('metadata_CH.txt') %>% as.data.frame()
stat_CH <- df_CH %>% group_by(dataset) %>% summarise(nMolecule = sum(UQ)) %>% as.data.frame()
colnames(stat_CH) <- c('dataset','nMolecule', dataset = 'CH-ATAC')

df_sci <- fread('metadata_sci-ATAC-seq.txt') %>% as.data.frame()
stat_sci <- data.frame(nMolecule = sum(df_sci$UQ), dataset = 'sci-ATAC-seq')

df_sci3 <- fread('metadata_sci3.txt') %>% as.data.frame()
stat_sci3 <- data.frame(nMolecule = sum(df_sci3$UQ), dataset = 'sci-ATAC-seq3')

df_CAT <- fread('metadata_CATlas.txt') %>% as.data.frame()
stat_CAT <- data.frame(nMolecule = sum(df_CAT$UQ), dataset = 'modified_sci-ATAC-seq')
  
### Tabula
metaDF_TH <- read.csv('Tabula_Sapiens_metadata.csv',row.names = 1)
TH_10X <- metaDF_TH[metaDF_TH$method == '10X',]
TH_10X_use <- TH_10X[,c(1,5)]
colnames(TH_10X_use) <- c('tissue','UQ')
TH_10X_use$tech <- 'Tabula Sapiens'
stat_human <- TH_10X_use %>% group_by(tech) %>% summarise(nMolecule = sum(UQ)) %>% as.data.frame()
colnames(stat_human) <- c('dataset','nMolecule')

metaDF_TM <- fread('TabulaMuris_metadata.txt') %>% as.data.frame()
TM_10X <- metaDF_TM[metaDF_TM$tech == 'Droplets',]
TM_10X_use <- TM_10X[,c(2,3)] 
colnames(TM_10X_use) <- c('tissue','UQ')
TM_10X_use$tech <- 'Tabula Muris'
stat_mouse <- TM_10X_use %>% group_by(tech) %>% summarise(nMolecule = sum(UQ)) %>% as.data.frame()
colnames(stat_mouse) <- c('dataset','nMolecule')

stat_all <- rbind(stat_UU, stat_CH, stat_sci, stat_sci3, stat_CAT, stat_mouse, stat_human)
stat_all
stat_all <- stat_all[order(stat_all$nMolecule, decreasing = F),]


stat_all$nMolecule <- as.numeric(stat_all$nMolecule)
stat_all$nMolecule_use <- stat_all$nMolecule/10^9
p <- ggplot(stat_all, aes(dataset, nMolecule_use, fill = dataset)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_text(aes(label = nMolecule), position=position_dodge(width = 0.9), size = 4, vjust = -0.25) +
  scale_fill_manual(values=c("#FFC700", "#A6CEE3", "#FB9FF0", "#30ABD6", "#FFF455", '#2A7FB0', "#E73429"))+ # "#F784EF"
  xlab("") + ylab("Total molecule number (x10^9)") + labs(title = "") +
  theme_classic() +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(size = 19, hjust = 0.5),
         axis.title = element_text(size = 18, color = "black"),
         axis.text = element_text(size = 18, color = "black"), axis.text.x = element_text(angle = 45, hjust = 1),
         legend.text = element_text(size=18))
print(p)
pdf(file = 'benchmark_scATAC&scRNA.pdf', w=11, h=6)
print(p)
dev.off()
