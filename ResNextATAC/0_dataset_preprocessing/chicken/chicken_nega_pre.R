##### set working dir ######
parser <- ArgumentParser()
parser$add_argument("-w", "--workdir", required=TRUE, help="absolute path of working directory")
args <- parser$parse_args()
workdir <- args$workdir
setwd(workdir)


###### load dataset ######
library(ArchR)
addArchRThreads(23)
load('./Gal6_archR.rda')
proj=readRDS('./proj_sample5w.rds')

###### obtain cell-peak binary matrix ######
peak_matrix <- getMatrixFromProject(proj, useMatrix = 'PeakMatrix', binarize = T)
pmat <- peak_matrix@assays@data@listData$PeakMatrix 
dim(pmat)
peak_info <- as.data.frame(proj@peakSet,row.names = 1:length(proj@peakSet))
peak_info$feature <- paste(peak_info$seqnames,peak_info$start,peak_info$end,sep = '-')

rownames(pmat) <- peak_info$feature
pmat[1:5,1:5]



###### filter sex chromosomes  ######
features_filter_chrZW <- peak_info$feature[!peak_info$seqnames %in% c('chrZ','chrW')]
pmat2 <- pmat[features_filter_chrZW,]



###### filter peaks accessible in less than 1% cells  ######
pmat3 <- pmat2[rowSums(pmat2)>dim(pmat2)[2]*0.01,]
dim(pmat3) 


###### fit pmat into h5ad format ######
library(SeuratDisk)
library(Seurat)
pmat3=pmat3[,meta_sample_use$barcode]
pmat3[1:10,1:10]
dim(pmat3)

seob <- CreateSeuratObject(counts = pmat3)
SaveH5Seurat(seob,filename = './tmp.h5Seurat',overwrite = T)
Convert('./tmp.h5Seurat',dest = './Chicken_5wCells_16wPeaks.h5ad',overwrite = T)
file.remove('./tmp.h5Seurat')


###### select negative samples ######
chromSizes <- Gal6_GenomeAnnotation$chromSizes
windows <- unlist(tile(chromSizes, width = 501)) 
windows


library(GenomicRanges)
peak=peak_matrix@rowRanges
peak
idx=queryHits(findOverlaps(windows,peak))
windows_negative=windows[-idx]
windows_negative

df=as.data.frame(windows_negative@ranges)
windows_negative=data.frame('chr'=windows_negative@seqnames,'start'=df$start,'end'=df$end)
windows_negative$width=windows_negative$end-windows_negative$start+1
head(windows_negative)
windows_negative=windows_negative[windows_negative$width==501,]
dim(windows_negative)
table(windows_negative$chr)
windows_negative=windows_negative[!windows_negative$chr%in%c('chrW','chrZ'),]
dim(windows_negative)
windows_negative_subset=windows_negative[sample(515806,80511),]
windows_negative=windows_negative_subset
windows_negative=windows_negative[!is.na(windows_negative$chr),]
table(windows_negative$chr)
windows_negative_peak=paste(windows_negative$chr,windows_negative$start,windows_negative$end,sep = '-')
fwrite(as.data.frame(windows_negative_peak),file = './nega_chicken.csv')
