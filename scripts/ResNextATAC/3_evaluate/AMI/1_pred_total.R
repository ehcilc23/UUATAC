###### set working dir ######
parser <- ArgumentParser()
parser$add_argument("-w", "--workdir", required=TRUE, help="absolute path of working directory")
args <- parser$parse_args()
workdir <- args$workdir
setwd(workdir)

mat=Matrix::readMM('pmat_pred_imputed_total.mtx')
dim(mat)
library(data.table)

peak=read.csv('./peak_total.csv')
head(peak)
peak$X0=gsub("b'",'',peak$X0)
peak$X0=gsub("'",'',peak$X0)
peak$X1=gsub("b'",'',peak$X1)
peak$X1=gsub("'",'',peak$X1)
peak$X2=gsub("b'",'',peak$X2)
peak$X2=gsub("'",'',peak$X2)
peak$X3=gsub("b'",'',peak$X3)
peak$X3=gsub("'",'',peak$X3)

peak$feature=paste(peak$X0,peak$X1,peak$X2,sep = '-')

rownames(mat)=peak$feature
dim(mat)
mat[1:5,1:5]
mat=as(mat,'dgCMatrix')


###### metadata pre ######
meta=fread('Female_mus_anno_5w.txt')
meta=as.data.frame(meta)
rownames(meta)=meta$barcode
colnames(mat)=meta$barcode
mat[1:5,1:5]

library(Matrix)
summary(rowSums(mat))
summary(colSums(mat))
dim(mat)


###### select peaks ######
idy=which(rowSums(mat)>50)
idy
mat1=mat[idy,]
dim(mat1)
#mat1=mat
#idx=which(rowSums(mat1)>100)


dim(mat1)


meta1=meta[colnames(mat1),]
mat1[1:5,1:5]
library(Signac)
library(Seurat)

###### signac pipeline ######
chrom_assay <- CreateChromatinAssay(
  counts = mat1,
  sep = c("-", "-")
)


obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = meta1,
  names.field = 1)
table(obj$anno)

obj=RunTFIDF(obj)
obj=FindTopFeatures(obj,min.cutoff = 'q0')
obj=RunSVD(obj)
DepthCor(obj,n = 50)

obj=FindNeighbors(obj,reduction = 'lsi',dims = 3:15)
obj=FindClusters(obj,resolution = 0.2,algorithm = 3,verbose = F)
table(obj$seurat_clusters)


obj=RunTSNE(obj,reduction = 'lsi',dims = 3:15,check_duplicates=F)

library(RColorBrewer)
library(ggsci)

DimPlot(obj,reduction = 'tsne',label = T,pt.size = 0.05)+NoLegend()
DimPlot(obj,group.by = 'lineage',label = T,reduction = 'tsne',cols = colorRampPalette(pal_npg(alpha = 0.6)(10))(10),pt.size = 0.05)

###### save ######
saveRDS(obj,file = 'obj_pred_total.rds')


tsne=as.data.frame(Embeddings(obj,reduction = 'tsne'))
head(tsne)
library(dplyr)

write.csv(tsne,'tsne_pred_total.csv')

write.csv(as.data.frame(obj@meta.data),'meta_cluster_pred_total.csv')

meta=fread('Female_mus_anno_5w.txt',header = T)
meta=as.data.frame(meta)
rownames(meta)=meta$barcode
meta1=meta[colnames(obj),]
write.csv(meta1,'meta_pred_total.csv')

