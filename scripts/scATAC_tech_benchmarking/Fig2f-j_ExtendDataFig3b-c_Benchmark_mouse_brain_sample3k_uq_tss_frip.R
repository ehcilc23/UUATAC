library(ggplot2)
library(Seurat)
library(dplyr)
library(plyr)
library(patchwork)

###### 10x-v1
load("proj_10x-v11_sample3k.rda")
proj_10xv11 = proj1_use
proj_10xv11$Platform = '10x-v1.1'

###### 10x-v2
load("proj_10x-v2_sample3k.rda")
proj_10xv2 = proj1_use
proj_10xv2$Platform = '10x-v2'

###### dscATAC
load("proj_dscATAC_sample3k.rda")
proj_dscATAC = proj1_use
proj_dscATAC$Platform = 'dscATAC'

###### HyDrop-ATAC
load("proj_HyDrop-ATAC_sample3k.rda")
proj_HyDropATAC = proj1_use
proj_HyDropATAC$Platform = 'HyDrop-ATAC'

###### Paired-seq
load("proj_Paired-seq_sample3k.rda")
proj_Pairedseq = proj1_use
proj_Pairedseq$Platform = 'Paired-seq'

###### sci-ATAC-seq
load("proj_sci-ATAC-seq_sample3k.rda")
proj_sciATACseq = proj1_use
proj_sciATACseq$Platform = 'sci-ATAC-seq'

###### SHARE-seq
load("proj_SHARE-seq_sample3k.rda")
proj_SHAREseq = proj1_use
proj_SHAREseq$Platform = 'SHARE-seq'

###### snATAC
load("proj_snATAC_sample3k.rda")
proj_snATAC = proj1_use
proj_snATAC$Platform = 'snATAC'

###### UUATAC
load("proj_UUATAC_sample3k.rda")
proj_UUATAC = proj1_use
proj_UUATAC$Platform = 'UUATAC'

df = as.data.frame(matrix(nrow=0,ncol=4)) 
for (i in c(proj_10xv11, proj_10xv2, proj_dscATAC, proj_HyDropATAC, proj_Pairedseq, proj_sciATACseq, proj_SHAREseq, proj_snATAC, proj_UUATAC))  
{
  df_temp <- data.frame( 'Sample' = i$Sample, 'TSS' = i$TSSEnrichment, 'UQ' = i$nFrags, 'FRIP' = i$FRIP, 'tech' = i$Platform)
  df <- rbind(df, df_temp)
}
head(df); dim(df); unique(df$tech)

color <- c("#9ED57B", "#5AB348", "#A67B5B", "#F4A261", "#9A9CC9", "#36BA98", "#FB9FF0", "#F3CA52", "#A6CEE3", "#DA7297", "#686D76", "#E73429", '#2A7FB0') 
names(color) <- c("10x-v1.1", "10x-v2", "dscATAC", "HyDrop-ATAC", "Paired-seq", "s3-ATAC", "CHATAC", "sciATACv3", "sci-ATAC-seq", "SHARE-seq", "snATAC", "UUATAC", "sci-ATAC-seq3")

################# UQ
df_UQ <- df %>% group_by(tech) %>% dplyr::summarize(median = median(UQ, na.rm=TRUE))
df_UQ <- dplyr::arrange(df_UQ, df_UQ[, c('median')])
df_UQ$tech

p_UQ <- ggplot(df, aes(tech, log10(UQ), fill = tech)) +
  geom_boxplot() +
  scale_fill_manual(values = color)+
  xlab("") + ylab("Log10 (Unique Fragments) in 3000 nuclei") + labs(title = "") +
  theme_bw() +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(size = 19, hjust = 0.5),
         axis.title = element_text(size = 18, color = "black"),
         axis.text = element_text(size = 18, color = "black"), axis.text.x = element_text(angle = 45, hjust = 1),
         legend.text = element_text(size=18)) +
  NoLegend()
print(p_UQ)

################ FRIP
df_FRIP <- df %>% group_by(tech) %>% dplyr::summarize(median = median(FRIP, na.rm=TRUE))
df_FRIP <- dplyr::arrange(df_FRIP, df_FRIP[, c('median')])
df_FRIP$tech

p_FRIP <- ggplot(df, aes(tech, FRIP, fill = tech)) +
  geom_boxplot() +
  scale_fill_manual(values=color)+
  xlab("") + ylab("FRiP in 3000 nuclei") + labs(title = "") +
  theme_bw() +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(size = 19, hjust = 0.5),
         axis.title = element_text(size = 18, color = "black"),
         axis.text = element_text(size = 18, color = "black"), axis.text.x = element_text(angle = 45, hjust = 1),
         legend.text = element_text(size=18)) +
  NoLegend()
print(p_FRIP)

################ Peaks
df_peak = as.data.frame(matrix(nrow=0, ncol=3)) 
for (i in c(proj_10xv11, proj_10xv2, proj_dscATAC, proj_HyDropATAC, proj_Pairedseq, 
            # proj_sciATACv3, 
            proj_sciATACseq, proj_SHAREseq, proj_snATAC, proj_UUATAC))  {
  df_temp <- data.frame( 'Sample' = unique(i$Sample), 'peaks' = length(i@peakSet@elementMetadata@listData$score))
  df_peak <- rbind(df_peak, df_temp)
}
df_peak$tech <- sub('_sample3k', '', df_peak$Sample)
df_peak$tech <- sub('10x-v11', '10x-v1.1', df_peak$tech)
head(df_peak); dim(df_peak)

p_peak <- ggplot(df_peak, aes(tech, peaks, fill = tech)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_text(aes(label = peaks), position=position_dodge(width = 0.9), size = 8, vjust = -0.25) + 
  scale_fill_manual(values=color)+
  xlab("") + ylab("Total Peaks in 3000 nuclei (x10^4)") + labs(title = "") +
  theme_classic() +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(size = 19, hjust = 0.5),
         axis.title = element_text(size = 18, color = "black"),
         axis.text = element_text(size = 18, color = "black"), axis.text.x = element_text(angle = 45, hjust = 1),
         legend.text = element_text(size=18)) +
  # ylim(7,30)+
  NoLegend()
print(p_peak)

pdf(file = '9tech_UQ_TSS_FRIP.pdf', w=18, h=12)
grid.arrange(p_UQ, p_FRIP, p_peak, nrow = 2)
dev.off()


################ Reads in Peaks
df_readspeak = as.data.frame(matrix(nrow=0,ncol=5)) 
for (i in c(proj_10xv11, proj_10xv2, proj_dscATAC, proj_HyDropATAC, proj_Pairedseq, proj_sciATACseq, proj_SHAREseq, proj_snATAC, proj_UUATAC))
{
  df_temp <- data.frame( 'Sample' = i$Sample, 'ReadsInPeaks' = i$ReadsInPeaks, 'ReadsInPromoter' = i$ReadsInPromoter, 'PromoterRatio' = i$PromoterRatio, 'tech' = i$Platform)
  df_readspeak <- rbind(df_readspeak, df_temp)
}
head(df_readspeak); dim(df_readspeak)

################# ReadsInPeaks
df_readspeak_peaks <- df_readspeak %>% group_by(tech) %>% dplyr::summarize(median = median(ReadsInPeaks, na.rm=TRUE))
df_readspeak_peaks <- dplyr::arrange(df_readspeak_peaks, df_readspeak_peaks[, c('median')])
df_readspeak_peaks$tech

p_readspeak_peaks <- ggplot(df_readspeak, aes(tech, log10(ReadsInPeaks), fill = tech)) +
  geom_boxplot() +
  scale_fill_manual(values = color)+
  xlab("") + ylab("Log10 (ReadsInPeaks) in 3000 nuclei") + labs(title = "") +
  theme_bw() +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(size = 19, hjust = 0.5),
         axis.title = element_text(size = 18, color = "black"),
         axis.text = element_text(size = 18, color = "black"), axis.text.x = element_text(angle = 45, hjust = 1),
         legend.text = element_text(size=18)) +
  NoLegend()
print(p_readspeak_peaks)


################# ReadsInPromoter
df_readspeak_promotor <- df_readspeak %>% group_by(tech) %>% dplyr::summarize(median = median(ReadsInPromoter, na.rm=TRUE))
df_readspeak_promotor <- dplyr::arrange(df_readspeak_promotor, df_readspeak_promotor[, c('median')])
df_readspeak_promotor$tech

p_readspeak_promotor <- ggplot(df_readspeak, aes(tech, log10(ReadsInPromoter), fill = tech)) +
  geom_boxplot() +
  scale_fill_manual(values = color)+
  xlab("") + ylab("Log10 (ReadsInPromoter) in 3000 nuclei") + labs(title = "") +
  theme_bw() +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(size = 19, hjust = 0.5),
         axis.title = element_text(size = 18, color = "black"),
         axis.text = element_text(size = 18, color = "black"), axis.text.x = element_text(angle = 45, hjust = 1),
         legend.text = element_text(size=18)) +
  NoLegend()
print(p_readspeak_promotor)


################# PromoterRatio
df_readspeak_ratio <- df_readspeak %>% group_by(tech) %>% dplyr::summarize(median = median(PromoterRatio, na.rm=TRUE))
df_readspeak_ratio <- dplyr::arrange(df_readspeak_ratio, df_readspeak_ratio[, c('median')])
df_readspeak_ratio$tech

p_readspeak_ratio <- ggplot(df_readspeak, aes(tech, PromoterRatio, fill = tech)) +
  geom_boxplot() +
  scale_fill_manual(values = color)+
  xlab("") + ylab("PromoterRatio in 3000 nuclei") + labs(title = "") +
  theme_bw() +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(size = 19, hjust = 0.5),
         axis.title = element_text(size = 18, color = "black"),
         axis.text = element_text(size = 18, color = "black"), axis.text.x = element_text(angle = 45, hjust = 1),
         legend.text = element_text(size=18)) +
  NoLegend()
print(p_readspeak_ratio)

pdf(file = '9tech_peaks_information.pdf', w=18, h=12)
grid.arrange(p_readspeak_peaks, p_readspeak_promotor, p_readspeak_ratio, nrow = 2)
dev.off()


#####################################################################################
################# TSS & UQ & FRIP lineplot
df = as.data.frame(matrix(nrow=0,ncol=4)) 
for (i in c(proj_10xv11, proj_10xv2, proj_dscATAC, proj_HyDropATAC, proj_Pairedseq, proj_sciATACseq, proj_SHAREseq, proj_snATAC, proj_UUATAC))  
{
  df_temp <- data.frame( 'Sample' = i$Sample, 'TSS' = i$TSSEnrichment, 'UQ' = i$nFrags, 'FRIP' = i$FRIP, 'tech' = i$Platform)
  df <- rbind(df, df_temp)
}
head(df); dim(df)

################# UQ & FRIP
p_UQ_FRIP <- ggplot(df, aes(log10(UQ), FRIP, color = tech)) +
  geom_smooth(method = 'lm',se=T,formula = y~x,fill = 'gray90') +
  scale_color_manual(values = color) +
  xlab("Log10 (Unique Fragments) in 3000 nuclei") + ylab("FRiP in 3000 nuclei") + labs(title = "") +
  theme_bw() +
  theme( plot.title = element_text(size = 19, hjust = 0.5),
         axis.title = element_text(size = 18, color = "black"),
         axis.text = element_text(size = 18, color = "black"), 
         legend.text = element_text(size=18)) +
  xlim(3.5, 6)
print(p_UQ_FRIP)

################# TSS & FRIP
p_TSS_FRIP <- ggplot(df, aes(TSS, FRIP, color = tech)) +
  geom_smooth(method = 'lm',se=T,formula = y~x,fill = 'gray90') +
  scale_color_manual(values = color) +
  xlab("TSS Enrichment in 3000 nuclei") + ylab("FRiP in 3000 nuclei") + labs(title = "") +
  theme_bw() +
  theme( plot.title = element_text(size = 19, hjust = 0.5),
    axis.title = element_text(size = 18, color = "black"),
    axis.text = element_text(size = 18, color = "black"), 
    legend.text = element_text(size=18))
print(p_TSS_FRIP)

pdf(file = '9tech_UQ_TSS_FRIP_lineplot.pdf', w=15.8, h=6)
grid.arrange(p_UQ_FRIP, p_TSS_FRIP, ncol = 2)
dev.off()


#####################################################################################
################# doublet score
ArrowFiles <- c(
  "10x-v11_sample3k.arrow",
  "10x-v2_sample3k.arrow",
  "dscATAC_sample3k.arrow",
  "HyDrop-ATAC_sample3k.arrow",
  "Paired-seq_sample3k.arrow",
 
  "sci-ATAC-seq_sample3k.arrow",
  "SHARE-seq_sample3k.arrow",
  "snATAC_sample3k.arrow",
  "UUATAC_sample3k.arrow"
)

df = as.data.frame(matrix(nrow=0, ncol=5)) 
for (i in 1:length(ArrowFiles)){
  proj1 <- ArchRProject(ArrowFiles = ArrowFiles[i], outputDirectory = "Output", copyArrows = FALSE)
  proj1 <- filterDoublets(proj1)
  df_temp <- data.frame( 'Sample' = proj1$Sample, 'TSS' = proj1$TSSEnrichment, 'UQ' = proj1$nFrags, 'DoubletEnrichment' = proj1$DoubletEnrichment)
  df <- rbind(df, df_temp)
}
head(df); dim(df)
df$tech <- sub('_sample3k', '', df$Sample)
df$tech <- sub('10x-v11', '10x-v1.1', df$tech)
unique(df$tech)

df_score <- df %>% group_by(tech) %>% dplyr::summarize(median_UQ = median(UQ, na.rm=TRUE),  
                                                         median_DoubletEnrichment = median(DoubletEnrichment, na.rm=TRUE),
                                                         mean_UQ = mean(UQ, na.rm=TRUE),
                                                         mean_DoubletEnrichment = mean(DoubletEnrichment, na.rm=TRUE))
head(df_score); dim(df_score)

p_UQ_DoubletEnrichment_mean <- ggscatter(df_score, 'mean_UQ', y = 'mean_DoubletEnrichment', color = 'tech', palette = color,
                                           label = 'tech', font.label = c(18, "plain"), size = 5) +
  xlab("Mean Unique Fragments") + ylab("Mean Doublet Enrichment") + labs(title = "") +
  theme_classic() +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(size = 19, hjust = 0.5),
         axis.title = element_text(size = 18, color = "black"),
         axis.text = element_text(size = 18, color = "black"),
         legend.text = element_text(size=18)) 
p_UQ_DoubletEnrichment_mean

pdf(file = '9tech_doublet.pdf', w=9, h=6)
print(p_UQ_DoubletEnrichment_mean)
dev.off()
