library(plyr)
library(ArchR)
library(ggplot2)
library(parallel)
library(patchwork)
library(ggpubr)
library(dplyr)
library(future)
plan()
plan("multisession", workers = 16)
plan()
addArchRThreads(threads = 16)
addArchRGenome("mm10")

######################## tech
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

color <- c("#9ED57B", "#5AB348", "#A67B5B", "#F4A261", "#9A9CC9", "#36BA98", "#FB9FF0", "#F3CA52", "#A6CEE3", "#DA7297", "#686D76", "#E73429", '#2A7FB0') 
names(color) <- c("10x-v1.1", "10x-v2", "dscATAC", "HyDrop-ATAC", "Paired-seq", "s3-ATAC", "CHATAC", "sciATACv3", "sci-ATAC-seq", "SHARE-seq", "snATAC", "UUATAC", "sci-ATAC-seq3")


######################## 
cell_use = as.data.frame(matrix(nrow=0, ncol=2)) 
for (i in 1:length(ArrowFiles)){
  proj1 <- ArchRProject(ArrowFiles = ArrowFiles[i], outputDirectory = "Output", copyArrows = FALSE)
  proj1 <- filterDoublets(proj1)
}
head(cell_use); dim(cell_use)
table(cell_use$Sample)  

######################## 
proj1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "Output",
  copyArrows = FALSE, #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
proj1
proj_use <- proj1[cell_use$cellname,]
proj_use
proj_use$Platform <- sub('_sample3k', '', proj_use$Sample)
proj_use$Platform <- sub('10x-v11', '10x-v1.1', proj_use$Platform)
table(proj_use$Platform)

# saveRDS(proj_use, 'proj_all9Tech_sample3k_deDoublet.rds')

######################## clustering
proj_rmbatch <- addIterativeLSI(
  ArchRProj = proj_use,
  useMatrix = 'TileMatrix',
  name = 'IterativeLSI',
  iterations = 3,
  clusterParams = list(
    resolution = c(2),
    sampleCells = 5000,
    n.start = 10
  ),
  varFeatures = 250000,
  filterQuantile = 0.95,
  dimsToUse = 1:15,
  force = TRUE
)

proj_rmbatch <- addHarmony(
  ArchRProj = proj_rmbatch,
  reducedDims = 'IterativeLSI',
  name = 'Harmony',
  groupBy = 'Sample'
)

proj_rmbatch <- addClusters(
  input = proj_rmbatch,
  reducedDims = 'Harmony',
  method = 'Seurat',
  name = 'Clusters',
  resolution = 0.5,
  maxClusters = NULL,
  force = TRUE
)

proj_rmbatch <- addUMAP(
  ArchRProj = proj_rmbatch,
  reducedDims = 'Harmony',
  name = 'UMAP',
  nNeighbors = 200,
  minDist = 0.4,
  metric = 'cosine',
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = proj_rmbatch, colorBy = 'cellColData', name = 'Clusters', embedding = 'UMAP',size = 0.1)
p2 <- plotEmbedding(ArchRProj = proj_rmbatch, colorBy = 'cellColData', name = 'Platform', embedding = 'UMAP',size = 0.1, pal = color)
pdf('Lr3_vF250000_Ld15_Fq0.95_cr0.5_deDoublets_nNeighbor200_rmbatch_9tech.pdf',width = 10,height = 7,onefile = F)
ggAlignPlots(p1, p2, type = "h")
dev.off()
saveRDS(proj_rmbatch,file = 'Lr3_vF250000_Ld15_Fq0.95_cr0.5_deDoublets_nNeighbor200_rmbatch_9tech.rds')


proj_rmbatch$anno <- recode(proj_rmbatch$Clusters,
                        'C1'='Dopaminergic neuron',
                        'C2'='Neuron',
                        'C3'='GABAergic neuron',
                        'C4'='Grabule neuron',
                        'C6'='Cholinergic',
                        'C10'='Medium spiny neuron',
                        'C5'='Glutamatergic neuron (VGLUT2)',
                        'C7'='Glutamatergic neuron (VGLUT2)',
                        'C8'='Glutamatergic neuron (VGLUT1)',
                        'C9'='Glutamatergic neuron (VGLUT1)',
                        'C11'='Glutamatergic neuron (VGLUT1)',
                        'C12'='Glutamatergic neuron (VGLUT1)',
                        'C13'='Glutamatergic neuron (VGLUT1)',
                        'C14'='Microglia',
                        'C15'='Perityte',
                        'C16'='Perityte',
                        'C17'='Myelinated oligodendrocyte',
                        'C18'='OPC',
                        'C19'='Astrocyte',
)

library(RColorBrewer)
color_anno <- as.character(t(as.data.frame(ArchRPalettes[1])))[1:13]
names(color_anno) <- c('Dopaminergic neuron', 'Neuron', 'GABAergic neuron', 'Grabule neuron', 'Cholinergic', 
                       'Medium spiny neuron', 'Perityte', 'Microglia', 'Myelinated oligodendrocyte', 
                       'Glutamatergic neuron (VGLUT2)', 'Glutamatergic neuron (VGLUT1)', 'OPC', 'Astrocyte')
color_anno
p1 <- plotEmbedding(ArchRProj = proj_rmbatch, colorBy = 'cellColData', name = 'anno', embedding = 'UMAP',size = 0.05, pal = color_anno)
p2 <- plotEmbedding(ArchRProj = proj_rmbatch, colorBy = 'cellColData', name = 'Platform', embedding = 'UMAP',size = 0.05, pal = color)
pdf('Lr3_vF250000_Ld15_Fq0.95_cr0.5_deDoublets_nNeighbor200_rmbatch_9tech_anno.pdf', width = 10, height = 7, onefile = F)
ggAlignPlots(p1, p2, type = "h")
dev.off()
saveRDS(proj_rmbatch,file = 'Lr3_vF250000_Ld15_Fq0.95_cr0.5_deDoublets_nNeighbor200_rmbatch_9tech.rds')


######################## DA peak
proj = list()
for (i in unique(proj_rmbatch$Sample)){
  proj[[i]] = proj_rmbatch[proj_rmbatch$Sample == i,]
  proj[[i]]@projectMetadata@listData[["outputDirectory"]] = paste0("/media/ggj/ggjlab_188/PeijingZhang/UUATAC/Benchmark_new/benchmark_sample3k_20240701/arrow_da_peak/")
  proj[[i]] <- addGroupCoverages(ArchRProj = proj[[i]], force = TRUE,
                                 groupBy = "anno", threads = 16)
  proj[[i]] <- addReproduciblePeakSet(
    ArchRProj =  proj[[i]], force = TRUE,
    groupBy = "anno", 
    maxPeaks = 1500000,
    pathToMacs2 = pathToMacs2)
  proj[[i]] <- addPeakMatrix(proj[[i]])
}

markersPeaks_list = list()
for(i in unique(proj_rmbatch$Sample)){
  table(proj[[i]]$anno)
  if ( i == "10x-v11_sample3k" || i == "10x-v2_sample3k"){
    proj[[i]] = proj[[i]][proj[[i]]$anno != 'Medium spiny neuron',]
  }else if( i == "sci-ATAC-seq_sample3k" ){
    proj[[i]] = proj[[i]][proj[[i]]$anno != 'GABAergic neuron',]
  }
  markersPeaks_list[[i]] = getMarkerFeatures(
    ArchRProj = proj[[i]], 
    useMatrix = "PeakMatrix", 
    groupBy = "anno",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon")
}

markerList = list()
for(i in unique(proj_rmbatch$Sample)){
  markerList[[i]] = getMarkers(markersPeaks_list[[i]], cutOff = "Pval < 0.05")
}
markerList_number = list()
for(i in unique(proj_rmbatch$Sample)){
  markerList_number[[i]] = dim(do.call(rbind, markerList[[i]]@listData))[1]
}

df_dapeak = data.frame(matrix(0, nrow = 9))
df_dapeak$tech = c('10x-v1.1', '10x-v2', 'dscATAC', 'HyDrop-ATAC', 'Paired-seq', 'sci-ATAC-seq', 'SHARE-seq', 'snATAC', 'UUATAC')
df_dapeak$DA_peaks = c( 275279, 296557, 60962, 148866, 3060, 212131, 16865, 55221, 855649 )
df_dapeak$DA_peaks <- as.numeric(df_dapeak$DA_peaks)
df_dapeak = df_dapeak[,-1]
df_dapeak
head(df_dapeak)

df_dapeak$tech <- factor(df_dapeak$tech, levels = c("Paired-seq", "SHARE-seq", "snATAC", "dscATAC", "HyDrop-ATAC", "sci-ATAC-seq", "10x-v1.1", "10x-v2", "UUATAC" ), ordered = TRUE)
df_dapeak$peaks_use <- df_dapeak$DA_peaks/10000
p_dapeak <- ggplot(df_dapeak, aes(tech, peaks_use, fill = tech)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_text(aes(label = DA_peaks), position=position_dodge(width = 0.9), size = 8, vjust = -0.25) + 
  scale_fill_manual(values=color)+
  xlab("") + ylab("DA peaks in 3000 nuclei (x10^4)") + labs(title = "") +
  theme_classic() +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title = element_text(size = 19, hjust = 0.5),
         axis.title = element_text(size = 18, color = "black"),
         axis.text = element_text(size = 18, color = "black"), axis.text.x = element_text(angle = 45, hjust = 1),
         legend.text = element_text(size=18)) +
  NoLegend()
print(p_dapeak)

pdf(file = '9tech_DAPeak.pdf', w=9, h=6)
print(p_dapeak)
dev.off()



############################ 
featureDF = topFeatures
mat <- .getPartialMatrix(
  ArrowFiles = ArrowFiles,
  featureDF = featureDF,
  useMatrix = useMatrix,
  cellNames = cellNames,
  doSampleCells = FALSE,
  threads = threads,
  verbose = FALSE
)

############ bin mean
proj_rmbatch$platCell = paste0(proj_rmbatch$Platform, '|', proj_rmbatch$anno)
df_platCell = table(proj_rmbatch$platCell) %>% as.data.frame()
df_platCell[which(df_platCell$Freq < 10),]
mat_use = as.matrix(mat) %>% as.data.frame()
rownames(mat_use) = paste0('bin_', 1:dim(mat_use)[1])
mat_use = t(mat_use) %>% as.data.frame()
mat_use[1:5,1:5]; dim(mat_use)
mat_use$platCell = proj_rmbatch$platCell

merge = aggregate(mat_use[1:20000], list(mat_use$platCell), mean)
rownames(merge) = merge[,1]
merge = merge[,-1]
merge_use = merge
merge_use = t(merge_use)
merge_use[1:5,1:5]

############ heatmap
library(lsa)
library(stringr)
df = cosine(merge_use)
df_use = df
df_use[1:5, 1:5]
rownames(df)
df_name = as.data.frame(rownames(df))
label = c("10x-v1.1|Astrocyte", "10x-v2|Astrocyte", "dscATAC|Astrocyte", "HyDrop-ATAC|Astrocyte", "Paired-seq|Astrocyte", "sci-ATAC-seq|Astrocyte", "SHARE-seq|Astrocyte", "snATAC|Astrocyte", "UUATAC|Astrocyte", "10x-v1.1|Cholinergic", "10x-v2|Cholinergic", "dscATAC|Cholinergic", "HyDrop-ATAC|Cholinergic", "Paired-seq|Cholinergic", "sci-ATAC-seq|Cholinergic", "SHARE-seq|Cholinergic", "snATAC|Cholinergic", "UUATAC|Cholinergic", "sci-ATAC-seq|Dopaminergic neuron", "UUATAC|Dopaminergic neuron", "10x-v1.1|GABAergic neuron", "10x-v2|GABAergic neuron", "dscATAC|GABAergic neuron", "HyDrop-ATAC|GABAergic neuron", "Paired-seq|GABAergic neuron", "sci-ATAC-seq|GABAergic neuron", "SHARE-seq|GABAergic neuron", "snATAC|GABAergic neuron", "10x-v1.1|Glutamatergic neuron (VGLUT1)", "10x-v2|Glutamatergic neuron (VGLUT1)", "dscATAC|Glutamatergic neuron (VGLUT1)", "HyDrop-ATAC|Glutamatergic neuron (VGLUT1)", "Paired-seq|Glutamatergic neuron (VGLUT1)", "sci-ATAC-seq|Glutamatergic neuron (VGLUT1)", "SHARE-seq|Glutamatergic neuron (VGLUT1)", "snATAC|Glutamatergic neuron (VGLUT1)", "UUATAC|Glutamatergic neuron (VGLUT1)", "10x-v1.1|Glutamatergic neuron (VGLUT2)", "10x-v2|Glutamatergic neuron (VGLUT2)", "dscATAC|Glutamatergic neuron (VGLUT2)", "HyDrop-ATAC|Glutamatergic neuron (VGLUT2)", "Paired-seq|Glutamatergic neuron (VGLUT2)", "sci-ATAC-seq|Glutamatergic neuron (VGLUT2)", "SHARE-seq|Glutamatergic neuron (VGLUT2)", "snATAC|Glutamatergic neuron (VGLUT2)", "UUATAC|Glutamatergic neuron (VGLUT2)", "Paired-seq|Grabule neuron", "sci-ATAC-seq|Grabule neuron", "SHARE-seq|Grabule neuron", "snATAC|Grabule neuron", "UUATAC|Grabule neuron", "10x-v1.1|Medium spiny neuron", "10x-v2|Medium spiny neuron", "dscATAC|Medium spiny neuron", "HyDrop-ATAC|Medium spiny neuron", "Paired-seq|Medium spiny neuron", "sci-ATAC-seq|Medium spiny neuron", "SHARE-seq|Medium spiny neuron", "snATAC|Medium spiny neuron", "10x-v1.1|Microglia", "10x-v2|Microglia", "dscATAC|Microglia", "HyDrop-ATAC|Microglia", "Paired-seq|Microglia", "sci-ATAC-seq|Microglia", "SHARE-seq|Microglia", "snATAC|Microglia", "UUATAC|Microglia", "10x-v1.1|Myelinated oligodendrocyte", "10x-v2|Myelinated oligodendrocyte", "dscATAC|Myelinated oligodendrocyte", "HyDrop-ATAC|Myelinated oligodendrocyte", "Paired-seq|Myelinated oligodendrocyte", "sci-ATAC-seq|Myelinated oligodendrocyte", "SHARE-seq|Myelinated oligodendrocyte", "snATAC|Myelinated oligodendrocyte", "UUATAC|Myelinated oligodendrocyte", "sci-ATAC-seq|Neuron", "UUATAC|Neuron", "10x-v1.1|OPC", "10x-v2|OPC", "dscATAC|OPC", "HyDrop-ATAC|OPC", "Paired-seq|OPC", "sci-ATAC-seq|OPC", "SHARE-seq|OPC", "snATAC|OPC", "UUATAC|OPC", "10x-v1.1|Perityte", "10x-v2|Perityte", "dscATAC|Perityte", "HyDrop-ATAC|Perityte", "Paired-seq|Perityte", "sci-ATAC-seq|Perityte", "SHARE-seq|Perityte", "snATAC|Perityte", "UUATAC|Perityte")
df_use <- df_use[order(factor(rownames(df_use) ,levels = label)),]
df_use <- df_use[,order(factor(colnames(df_use) ,levels = label))]

library(pheatmap)
library(RColorBrewer)
p_heatmap <- pheatmap(df_use, scale = 'none', border = 'white',
                      legend = T,legend_breaks = c(0,0.2,0.4,0.6,0.8,1),
                      treeheight_row = 20,treeheight_col = 20,
                      cluster_cols = T,cluster_rows = T,
                      clustering_method = 'ward.D2',
                      color = colorRampPalette(rev(brewer.pal(n=7,name = 'RdBu')))(1000))
print(p_heatmap)
pdf(file = '9tech_celltype_heatmap.pdf',width = 8.5,height = 8)
print(p_heatmap)
dev.off()


############ tree
library(grDevices)
library(dendextend)

total.dist <- as.dist(1-df_use)
total.tree <- hclust(total.dist,method = 'ward.D')
total.tree <- as.dendrogram(total.tree)

celltype_color = c(rep("#D51F26", 9), rep("#272E6A", 9), rep("#208A42", 2), rep("#89288F", 8), rep("#F47D2B", 9), 
                    rep("#FEE500", 9), rep("#8A9FD1", 5), rep("#C06CAB", 8), rep("#E6C2DC", 9),rep("#90D5E4", 9),
                    rep("#89C75F", 2), rep("#F37B7D", 9), rep("#9983BD", 9) )
names(celltype_color) = rownames(df_use)

dend<-total.tree %>%
  color_branches(k=6) %>%
  set("labels_cex",0.6)%>%
  set("leaves_pch",19)%>%
  set("leaves_col",celltype_color) %>%
  set("labels_colors",celltype_color)

cluster.list <- unique(names(celltype_color))
Celltype_bar<-matrix(nrow=length(cluster.list),ncol=length(cluster.list))
for (i in 1:length(cluster.list)){
  Celltype_bar[i,]<-ifelse(names(celltype_color) == as.character(cluster.list[i]),as.character(celltype_color),"white")	
}
Celltype_bar<-t(Celltype_bar)
dend1<- color_branches(dend,k=6)
plot(dend)

pdf(file = '9tech_celltype_tree.pdf')
par(mar=c(7,3,7,3)+0.1,xpd=NA)
plot(dend)
dev.off()
