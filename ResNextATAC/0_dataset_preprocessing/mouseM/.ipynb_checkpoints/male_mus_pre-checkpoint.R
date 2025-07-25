##### set working dir ######
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

###### remove blacklist ######
library(ArchR)
ArchR::addArchRGenome('mm10')
blacklist=genomeAnnoMm10$blacklist
id=queryHits(findOverlaps(peak_matrix@rowRanges,blacklist))
id=unique(id)
peak_info=peak_info[-id,]

###### filter sex chromosomes  ######
features_filter_chrXY <- peak_info$feature[peak_info$seqnames != 'chrX']
pmat2 <- pmat[features_filter_chrXY,]



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
colnames(seob)[1:5]
SaveH5Seurat(seob,filename = './Mouse_male.h5Seurat',overwrite = T)
Convert('./Mouse_male.h5Seurat',dest = './Mouse_Male_5wCells_15wPeaks.h5ad',overwrite = T)
file.remove('./Mouse_male.h5Seurat')

###### select negative samples ######
library(BSgenome.Mmusculus.UCSC.mm10)
genome <- BSgenome.Mmusculus.UCSC.mm10
chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
windows <- unlist(tile(chromSizes, width = 501)) 
windows


library(GenomicRanges)
idx=queryHits(findOverlaps(windows,peak_matrix@rowRanges))
windows_negative=windows[-idx]
windows_negative
idy=queryHits(findOverlaps(windows_negative,blacklist))
windows_negative=windows_negative[-idy]
windows_negative

df=as.data.frame(windows_negative@ranges)
windows_negative=data.frame('chr'=windows_negative@seqnames,'start'=df$start,'end'=df$end)
windows_negative$width=windows_negative$end-windows_negative$start+1
head(windows_negative)
windows_negative=windows_negative[windows_negative$width==501,]
dim(windows_negative)
table(windows_negative$chr)
windows_negative=windows_negative[!windows_negative$chr%in%c('chrM','chrY','chrX'),]
dim(windows_negative)


windows_negative_subset=windows_negative[sample(2671109,79214),]

windows_negative_subset=windows_negative_subset[!is.na(windows_negative_subset$chr),]
table(windows_negative_subset$chr)
windows_negative_subset=paste(windows_negative_subset$chr,windows_negative_subset$start,windows_negative_subset$end,sep = '-')
fwrite(as.data.frame(windows_negative_subset),file = 'nega_malemus.csv')