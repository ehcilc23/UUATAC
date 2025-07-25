library(data.table)
library(dplyr)
library(reshape)
library(ArchR)

## UU
proj_UU_musM <- readRDS('proj_raw_UUATAC.rds')
proj_UU_musM <- proj_UU_musM[proj_UU_musM$TSSEnrichment>7&proj_UU_musM$nFrags>1000,]
df_UU <- getCellColData(ArchRProj = proj_UU_musM, select = c('nFrags','TSSEnrichment')) %>% as.data.frame()
df_UU$tech <- 'UUATAC'

## 10x-v1.1
proj1 <- readRDS('proj_raw_10x-v11.rds')
proj1 <- proj1[proj1$TSSEnrichment>7&proj1$nFrags>1000,]
df_10x_v11 <- getCellColData(ArchRProj = proj1, select = c('nFrags','TSSEnrichment')) %>% as.data.frame()
df_10x_v11$tech <- '10x-v1.1'

## 10x-v2
proj1 <- readRDS('proj_raw_10x-v2.rds')
proj1 <- proj1[proj1$TSSEnrichment>7&proj1$nFrags>1000,]
df_10x_v2 <- getCellColData(ArchRProj = proj1, select = c('nFrags','TSSEnrichment')) %>% as.data.frame()
df_10x_v2$tech <- '10x-v2'

## dsciATAC
proj1 <- readRDS('proj_raw_dsciATAC.rds')
df_dsci <- getCellColData(ArchRProj = proj1, select = c('nFrags','TSSEnrichment')) %>% as.data.frame()
df_dsci <- df_dsci[df_dsci$nFrags>1000 & df_dsci$TSSEnrichment>7,]
df_dsci$tech <- 'dsciATAC'

## HyDrop
proj1 <- readRDS('proj_raw_HyDrop-ATAC.rds')
proj1 <- proj1[proj1$TSSEnrichment>7&proj1$nFrags>1000,]
df_HyDrop <- getCellColData(ArchRProj = proj1, select = c('nFrags','TSSEnrichment')) %>% as.data.frame()
df_HyDrop$tech <- 'HyDrop-ATAC'

## Paired-seq
proj1 <- readRDS('proj_raw_Paired-seq.rds')
proj1 <- proj1[proj1$TSSEnrichment>7&proj1$nFrags>1000,]
df_paired <- getCellColData(ArchRProj = proj1, select = c('nFrags','TSSEnrichment')) %>% as.data.frame()
df_paired$tech <- 'Paired-seq'

## s3ATAC
proj1 <- readRDS('proj_raw_s3ATAC.rds')
proj1 <- proj1[proj1$TSSEnrichment>7&proj1$nFrags>1000,]
df_s3ATAC <- getCellColData(ArchRProj = proj1, select = c('nFrags','TSSEnrichment')) %>% as.data.frame()
df_s3ATAC$tech <- 's3-ATAC'

## SHARE-seq
proj1 <- readRDS('proj_raw_SHARE-seq.rds')
proj1 <- proj1[proj1$TSSEnrichment>7&proj1$nFrags>1000,]
df_SHARE <- getCellColData(ArchRProj = proj1, select = c('nFrags','TSSEnrichment')) %>% as.data.frame()
df_SHARE$tech <- 'SHARE-seq'

## snATAC
proj1 <- readRDS('proj_raw_snATAC.rds')
proj1 <- proj1[proj1$TSSEnrichment>7&proj1$nFrags>1000,]
df_snATAC <- getCellColData(ArchRProj = proj1, select = c('nFrags','TSSEnrichment')) %>% as.data.frame()
df_snATAC$tech <- 'snATAC'

## sci-ATAC-seq
df_sci <- fread('cellpass_QC_sci_batch_1.txt') %>% as.data.frame()
rownames(df_sci) <- df_sci$V1
df_sci <- df_sci[, c('nFrags', 'TSSEnrichment')]
df_sci <- df_sci[df_sci$nFrags>1000 & df_sci$TSSEnrichment>7,]
df_sci$tech <- 'sci-ATAC-seq'

## sci-ATAC-seq3
df_sci3 <- fread('cellpass_QC_sci3_humanFetal_batch_2.txt') %>% as.data.frame()
rownames(df_sci3) <- df_sci3$V1
df_sci3 <- df_sci3[, c('nFrags', 'TSSEnrichment')]
df_sci3 <- df_sci3[df_sci3$nFrags>1000 & df_sci3$TSSEnrichment>7,]
df_sci3$tech <- 'sci-ATAC-seq3'


## merge all techniques & plot --------------  
pal <- c("#9ED57B", "#5AB348", "#A67B5B", "#F4A261", "#9A9CC9", "#36BA98", "#FB9FF0",
         "#F3CA52", "#2A7FB0", "#A6CEE3", "#DA7297", "#686D76", "#E73429")
names(pal) <- c("10x-v1.1", "10x-v2", "dsciATAC", "HyDrop-ATAC", "Paired-seq", "s3-ATAC", "CHATAC",
                  "sciATACv3", "sci-ATAC-seq3", "sci-ATAC-seq", "SHARE-seq", "snATAC", "UUATAC")

df <- rbind(df_UU, df_10x_v11, df_10x_v2, df_dsci, df_HyDrop, df_paired, df_s3ATAC, df_SHARE, df_snATAC, df_sci, df_sci3 )
table(df$tech)
df$UQxTSS <- df$nFrags * df$TSSEnrichment

stat.df <- df %>% group_by(tech) %>% summarise(IC = sum(UQxTSS)) %>% as.data.frame()
stat.df <- stat.df[order(stat.df$IC,decreasing = F),]

p1 <- ggbarplot(stat.df,
                x = "tech", y = "IC",
                ylab = "Information content (IC)",
                fill = "tech", 
                color = "#000000",
                palette = pal,
                label = TRUE,
                position = position_dodge(0.8)) + 
  labs(title = 'Information content (accumulated TSSxUQ) comparisons of scATAC techniques')
p2 <- ggpar(p1, legend = 'right', x.text.angle = 90)
p2
pdf('benchmark_IC.pdf',width = 10,height = 7)
print(p2)
dev.off()
