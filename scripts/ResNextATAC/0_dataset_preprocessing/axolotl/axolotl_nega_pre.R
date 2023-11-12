###### set working dir ######
parser <- ArgumentParser()
parser$add_argument("-w", "--workdir", required=TRUE, help="absolute path of working directory")
args <- parser$parse_args()
workdir <- args$workdir
setwd(workdir)


###### load dataset ######
library(ArchR)
addArchRThreads(23)
proj=readRDS('./proj_sample5w.rds')

###### obtain cell-peak binary matrix ######
peak_matrix <- getMatrixFromProject(proj, useMatrix = 'PeakMatrix', binarize = T)
pmat <- peak_matrix@assays@data@listData$PeakMatrix 
dim(pmat)
peak_info <- as.data.frame(proj@peakSet,row.names = 1:length(proj@peakSet))
peak_info$feature <- paste(peak_info$seqnames,peak_info$start,peak_info$end,sep = '-')

rownames(pmat) <- peak_info$feature
pmat[1:5,1:5]

###### filter peaks accessible in less than 1% cells  ######
pmat2=pmat
pmat3 <- pmat2[rowSums(pmat2)>dim(pmat2)[2]*0.01,]
dim(pmat3) 


###### fit pmat into h5ad format ######
library(SeuratDisk)
library(Seurat)
pmat3=pmat3[,meta_sample_use$barcode]
pmat3[1:5,1:5]
dim(pmat3)

seob <- CreateSeuratObject(counts = pmat3)
SaveH5Seurat(seob,filename = './tmp.h5Seurat',overwrite = T)
Convert('./tmp.h5Seurat',dest = './Axolotl_5wCells_8wPeaks.h5ad',overwrite = T)
file.remove('./tmp.h5Seurat')

###### rename positive samples ######
postive_peak=fread('postive_peak.csv',header = T)
head(postive_peak)
postive_peak$peak=gsub('p-1-','p_1-',postive_peak$peak)
postive_peak$peak=gsub('p-2-','p_2-',postive_peak$peak)
postive_peak$peak=gsub('p-3-','p_3-',postive_peak$peak)
postive_peak$peak=gsub('q-1-','q_1-',postive_peak$peak)
postive_peak$peak=gsub('q-2-','q_2-',postive_peak$peak)
postive_peak$peak=gsub('q-3-','q_3-',postive_peak$peak)
postive_peak$peak=gsub('q-4-','q_4-',postive_peak$peak)
fwrite(postive_peak,file = 'postive_peak_new.csv')

###### select negative samples ######
chromSizes=fread('/media/ggj/Guo-4T-AI/FYT/UU/ref/AmexG_v6/ax.chrom.size')
chromSizes=GRanges(seqnames = chromSizes$V1,IRanges(1,chromSizes$V2))

chromSizes
windows <- unlist(tile(chromSizes, width = 501)) 
windows@seqnames


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
dim(windows_negative)

windows_negative_subset=windows_negative[sample(36733934,41318),]
windows_negative_subset=windows_negative_subset[!is.na(windows_negative_subset$chr),]
table(windows_negative_subset$chr)
windows_negative=windows_negative_subset
windows_negative_peak=paste(windows_negative$chr,windows_negative$start,windows_negative$end,sep = '-')
windows_negative_peak[1:5]
fwrite(as.data.frame(windows_negative_peak),file = './nega_ax.csv')
