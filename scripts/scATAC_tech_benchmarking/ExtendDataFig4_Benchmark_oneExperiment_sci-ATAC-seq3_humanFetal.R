library(data.table)
library(dplyr)
library(reshape)
library(ArchR)

############# load in proj
proj.ls <- list.files('./proj_raw_sci3humanFetal',pattern = '*.rds',full.names = T)
proj.ls

meta.df <- data.frame()
for (i in 1:length(proj.ls)) {
  message(i)
  proj_tmp <- readRDS(proj.ls[i])
  proj_tmp <- proj_tmp[proj_tmp$TSSEnrichment>7 & proj_tmp$nFrags>1000,]
  df <- getCellColData(ArchRProj = proj_tmp, select = c('nFrags','TSSEnrichment')) %>% as.data.frame()
  df$`log10(nFrags)` <- log10(df$nFrags)
  
  meta.df <- rbind(meta.df,df)
  rm(proj_tmp);gc()
}
write.table(meta.df,file = 'cellpass_QC_sci3_humanfetal.txt',sep = '\t',quote = F)

df <- fread('cellpass_QC_sci3_humanfetal.txt') %>% as.data.frame()
head(df); dim(df)
rownames(df) <- df$V1
df$sample <- sub('\\#.*', '', df$V1)
table(df$sample)

df_batch_1 <- df[grepl('sample_1_adrenal|sample_2_thymus|sample_3_kidney|sample_4_lung|sample_5_kidney|sample_6_brain|sample_7_liver|sample_8_lung|sample_9_liver|sample_10_muscle|sample_11_brain|sample_12_heart|sample_13_placenta|sample_14_heart|sample_15_placenta|sample_16_adrenal|sample_17_cerebellum|sample_18_eye|sample_19_intestine|sample_20_kidney|sample_21_intestine|sample_22_kidney|sample_23_muscle', df$sample) ,] 
head(df_batch_1); dim(df_batch_1)
df_batch_1$batch <- 'batch_1'
# write.table(df_batch_1, file = 'cellpass_QC_sci3_humanfetal_batch_1.txt', sep = '\t',quote = F, row.names = F)

df_batch_2 <- df[grepl('sample_24_standard|sample_25_kidney|sample_26_placenta|sample_27_adrenal|sample_28_intestine|sample_29_placenta|sample_30_placenta|sample_31_adrenal|sample_32_heart|sample_33_lung|sample_34_kidney|sample_35_liver|sample_36_brain|sample_37_liver|sample_38_lung|sample_39_heart|sample_40_liver|sample_41_thymus|sample_42_heart|sample_43_liver|sample_44_lung|sample_45_brain|sample_46_liver|sample_47_lung|sample_48_standard', df$sample) ,]
head(df_batch_2); dim(df_batch_2)
df_batch_2$batch <- 'batch_2'
# write.table(df_batch_2, file = 'cellpass_QC_sci3_humanfetal_batch_2.txt', sep = '\t',quote = F, row.names = F)

df_batch_3 <- df[grepl('sample_49_stomach|sample_50_bonemarrow|sample_51_gonad|sample_52_pancreas|sample_53_eye|sample_54_thymus|sample_55_eye|sample_56_spleen|sample_57_spleen|sample_58_cerebellum|sample_59_bonemarrow|sample_60_thymus|sample_61_pancreas|sample_62_stomach|sample_63_gonad|sample_64_brain|sample_65_kidney|sample_66_brain|sample_67_kidney|sample_68_lung|sample_69_brain|sample_70_lung|sample_71_brain|sample_72_standard', df$sample) ,]
head(df_batch_3); dim(df_batch_3)
df_batch_3$batch <- 'batch_3'
# write.table(df_batch_3, file = 'cellpass_QC_sci3_humanfetal_batch_3.txt', sep = '\t',quote = F, row.names = F)


############# select one experiment
df <- df_batch_2[,2:4]
head(df); dim(df)

UQ_cutoff = 1000, TSS_cutoff = 7
df_sub <- df[df$nFrags>UQ_cutoff & df$TSSEnrichment>TSS_cutoff,]
dim(df_sub) 

UQinfo <- as.numeric(summary(df_sub$`log10(nFrags)`))
UQinfo <- data.frame(info = c('Min','1st Quantile','Median','Mean','3rd Quantile','Max'),UQinfo)
MedianUQ = round(10^(UQinfo$UQinfo[3]))

TSSinfo <- as.numeric(summary(df_sub$TSSEnrichment))
TSSinfo <- data.frame(info = c('Min','1st Quantile','Median','Mean','3rd Quantile','Max'),TSSinfo)
MedianTSS = round(TSSinfo$TSSinfo[3],2)

tissue <- 'sci-ATAC-seq3'
p <- ggPoint(
  x = df[,3],
  y = df[,2],
  colorDensity = TRUE,
  rastr = T,
  pal = colorRampPalette(viridis(10,option = 'D')[1:10])(256),
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
) + geom_hline(yintercept = TSS_cutoff, lty = "dashed") +
  geom_vline(xintercept = log10(UQ_cutoff), lty = "dashed") +
  labs(title = paste0('scATAC-QC of ',tissue),
       tag = paste0('nCells Pass Filter = ',nrow(df_sub),'\n',
                    'Median Frags = ', MedianUQ,'\n',
                    'Median TSS Enrichment = ',MedianTSS,'\n',
                    'Cut-off: TSS>',TSS_cutoff,',UQ>',UQ_cutoff)) +
  theme(plot.title = element_text(size = 15,hjust = 0.5),
        plot.tag = element_text(size = 13))
p

pdf(file = paste0('Figure_oneExperiment_',tissue,'_TSS1_batch_3','.pdf'), width = 8, height = 8)
print(p)
dev.off()
