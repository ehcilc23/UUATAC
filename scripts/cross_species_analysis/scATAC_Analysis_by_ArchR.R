###### scATAC workflow of one species ######
library(argparse)
library(dplyr)
library(data.table)
library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Seurat)


parser <- ArgumentParser()
parser$add_argument("-w", "--workdir", required=TRUE, help="absolute path of working directory")
args <- parser$parse_args()
workdir <- args$workdir
setwd(workdir)
tissue <- basename(getwd()) # the name of working directory is the tissue name

## set dir
ProjDir <- file.path(getwd(),'Proj')
FragmentsDir <- file.path(getwd(),'Fragments')
ArrowDir <- file.path(getwd(),'ArrowFiles_bin500')

## set params
addArchRThreads(threads = 1)
addArchRGenome("mm10")



## ------------ 1 Create ArrowFiles and project ------------
inputFiles <- list.files(path = FragmentsDir, pattern = '*sort.bed.gz$',full.names = T)
names(inputFiles) <- tissue
inputFiles

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 0,
  minFrags = 1,
  maxFrags = 1000000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  TileMatParams = list(tileSize = 500),  # bin size
  # GeneScoreMatParams = list(tileSize = 2000),  # gene body upstream/downstream 2k
  threads = 1 # must be set to 1
)

doubScore <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

ArrowFiles <- list.files(path = ArrowDir, pattern = '*arrow$', full.names = T)
ArrowFiles

proj1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "Output",
  copyArrows = FALSE, #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
proj1
proj1$tissue <- tissue
saveRDS(proj1, file = file.path(ProjDir,'proj_raw.rds'))



## ------------ 2.1 Filter TSS/UQ ------------
df <- getCellColData(ArchRProj = proj1, select = c('nFrags','TSSEnrichment'))
df$`log10(nFrags)` <- log10(df$nFrags)
dim(df)
p <- ggPoint(
  x = df[,3],
  y = df[,2],
  colorDensity = TRUE,
  continuousSet = "sambaNight") +
  geom_hline(yintercept = 13, lty = "dashed") +
  geom_vline(xintercept = 4, lty = "dashed")
p


## set cut-off
UQ_cutoff = 10000
TSS_cutoff = 12

df_sub <- df[df$nFrags>UQ_cutoff & df$TSSEnrichment>TSS_cutoff,]
dim(df_sub)
UQinfo <- as.numeric(summary(df_sub$`log10(nFrags)`))
UQinfo <- data.frame(info = c('Min','1st Quantile','Median','Mean','3rd Quantile','Max'),UQinfo)
MedianUQ = round(10^(UQinfo$UQinfo[3]))

TSSinfo <- as.numeric(summary(df_sub$TSSEnrichment))
TSSinfo <- data.frame(info = c('Min','1st Quantile','Median','Mean','3rd Quantile','Max'),TSSinfo)
MedianTSS = round(TSSinfo$TSSinfo[3],2)

p <- ggPoint(
  x = df[,3],
  y = df[,2],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
) + geom_hline(yintercept = TSS_cutoff, lty = "dashed") +
  geom_vline(xintercept = log10(UQ_cutoff), lty = "dashed") +
  labs(title = paste0('UU scATAC-QC of ',tissue),
       tag = paste0('nCells Pass Filter = ',nrow(df_sub),'\n',
                    'Median Frags = ', MedianUQ,'\n',
                    'Median TSS Enrichment = ',MedianTSS,'\n',
                    'Cut-off: TSS>',TSS_cutoff,',UQ>',UQ_cutoff)) +
  theme(plot.title = element_text(size = 15,hjust = 0.5),
        plot.tag = element_text(size = 13))
p


QC_file_name = paste0('QC-TSS',TSS_cutoff,'_UQ',UQ_cutoff,'.pdf')
pdf(file = file.path('./Output', QC_file_name),width = 8,height = 8)
print(p)
dev.off()

## filter cells that passed QC
cellsPass <- rownames(df_sub)
proj1 <- proj1[cellsPass,]



## ------------ 2.2 Filter doublets ------------
proj2 <- addDoubletScores(input = proj1,
                          dimsToUse = 1:7,
                          knnMethod = 'LSI',
                          LSIParams = list(iterations = 1,
                                           varFeatures = 3000,
                                           dimsToUse = 1:7,
                                           filterQuantile = 0.95),
                          force = T,verbose = T,threads = 1)


proj2 <- filterDoublets(proj2)
proj2 <- filterDoublets(proj2,cutEnrich = 0)


# save 
Filtered_file_name <- paste0('proj_TSS',TSS_cutoff,'_UQ',UQ_cutoff,'_deDoublets.rds')
saveRDS(proj2,file = file.path(ProjDir,Filtered_file_name))



## ------------ 3 Clustering and Annotating ------------
## load the project of all tissues for one species
proj_raw <- readRDS('/file/path/prefix/proj_Global.rds')

## packed clutering procedure
archr_cluster = function(proj_raw, 
                         LSI = 2,LSI_res = 0.2,sampleCells = 10000,varFeatures = 25000,LSI_dims = 10,quantile = 0.995, 
                         cluster_res = 0.8, nOutlier = 5, maxClusters = NULL)
{
  project <- addIterativeLSI(
    ArchRProj = proj_raw,
    useMatrix = 'TileMatrix',
    name = 'IterativeLSI',
    iterations = LSI,
    clusterParams = list(
      resolution = c(LSI_res),
      sampleCells = sampleCells,
      n.start = 10
    ),
    varFeatures = varFeatures, # set varFeatures
    dimsToUse = 1:LSI_dims, # set dims to use
    filterQuantile = quantile, # filter features
    outlierQuantiles = c(0.05,0.95), # filter cells
    force = TRUE
  )
  
  project <- addClusters(
    input = project,
    reducedDims = 'IterativeLSI',
    method = 'Seurat',
    name = 'Clusters',
    resolution = cluster_res, # set res
    maxClusters = maxClusters, 
    nOutlier = nOutlier, # set the minimal number of cells
    force = TRUE
  )
  
  project <- addUMAP(
    ArchRProj = project,
    reducedDims = 'IterativeLSI',
    name = 'UMAP',
    nNeighbors = 40,
    minDist = 0.4,
    metric = 'cosine',
    force = TRUE
  )
  return(project)
}


## run in loop, to optimize UMAP visualization
picDir <- '/file/path/prefix'

varFeatures.list = c(25000)
LSI_dims.list = c(20,25,30)
quantile.list = c(0.995,0.99)

for (varFeatures in varFeatures.list) {
  for (LSI_dims in LSI_dims.list) {
    for (quantile in quantile.list) {
      # set file name
      pic_name = paste0('Lr',2,'_','vF',varFeatures,'_','Ld',LSI_dims,'_','Fq',quantile,'_','cr',0.8,'_deDoublets','.pdf')
      proj_name = paste0('Lr',2,'_','vF',varFeatures,'_','Ld',LSI_dims,'_','Fq',quantile,'_','cr',0.8,'_deDoublets','.rds')
      
      # demensionality reduction & UMAP
      proj <- archr_cluster(proj_raw = proj_raw, LSI_res = 2, varFeatures = varFeatures, LSI_dims = LSI_dims, quantile = quantile, cluster_res = 0.8)
      
      # save plot
      p1 <- plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'Clusters', embedding = 'UMAP')
      p2 <- plotEmbedding(ArchRProj = proj, colorBy = 'cellColData', name = 'Sample', embedding = 'UMAP')
      pdf(file = file.path(picDir,pic_name),width = 10,height = 7)
      ggAlignPlots(p1, p2, type = "h")
      dev.off()
      
      # save proj
      saveRDS(proj,file = file.path(ProjDir,proj_name))
    }
  }
}



## ------------ 4 Get markerlist and Assign celltype annotation ------------
markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = 'GeneScoreMatrix',
  groupBy = 'Clusters',
  bias = c('TSSEnrichment','log10(nFrags)'),
  testMethod = 'wilcoxon'
  # useGroups = 'C11', #
  # bgdGroups = 'C7' #
)
markersGS
markerList <- getMarkers(markersGS, cutOff = 'FDR <= 0.01 & Log2FC >= 1')
write.csv(markerList, file = './Global_markerList.csv',quote = F) # export and check the top marker genes to assign annotation

markerGenes <- c('Cd8b1','Cd3g',# thymocyte
                 'Il7r', # T cell
                 'Cd79b','Cd79a',# B cell
                 'Vpreb3', # preB
                 'Jchain','Xbp1', # plasma
                 'C1qa','Csf1r', # macro
                 'Marco','Mcemp1',# Alveolar macrophage
                 'S100a8', # myeloid
                 'Hbb-y','Hbb-bs', # erythroid
                 'Fat2','Zic2', # granule
                 'Slc32a1','Gad2',# GABA
                 'Lin28b', # Lin28b
                 'Th','Phox2b', # chromaffin
                 'Sox10','Plp1', # oligo
                 'Gfap','Agt',# Astrocyte
                 'Prl', # Lactotrope
                 'Pomc','Fshb', # Corticotrope/Gonadotrope
                 'Gh', # Somatotrope
                 'Six1','Eya1', # pitui stem
                 'Acta1','Tnni2', # skeletal
                 'Tnnt2','Tnni3', # cardiac
                 'Myh11','Acta2', # smc
                 'Col1a1','Col6a1', # fibro
                 'Adipoq','Adrb3',# Adipo
                 'Kdr','Flt1',# endo
                 'Foxl2','Prss35',# granulosa
                 'Nr5a1','Star',# adrenal cortex
                 'Chga','Isl1', # endocrine , + CHGA
                 'Nkx6-1','Hhex', # pancreas secretory
                 'Cpa1','Cpa2', # Acinar
                 'Tff1','Gkn1', # foveolar
                 'Muc6', # Mucous neck cell
                 'Pgc','Gif', # chief
                 'Atp4b','Atp4a', # parietal
                 'Cdx1','Cdx2',# enterocyte
                 'Slc22a12','Slc6a18', # PT
                 'Slc12a1', 'Atp6v1b1', #dist and collect
                 'Fga','Alb', # hepa  
                 'Sftpb','Scnn1a', # pulm epi
                 'Upk1b','Upk1a', # Urothelial cell
                 'Trp63','Krt14',# Basal
                 'Trpv6','Padi1') # uterus

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = 'FDR <= 0.01 & Log2FC >= 1',
  labelMarkers = markerGenes,
  transpose = TRUE,
  returnMatrix = FALSE
)
pdf(file = file.path(picDir,'MarkersHeatmap.pdf'), width = 10, height = 7)
ComplexHeatmap::draw(heatmapGS,
                     heatmap_legend_side = 'bot', annotation_legend_side = 'bot',
                     # row_order = reorder_row,
                     cluster_rows = TRUE)
dev.off()

## MAGIC(imputation)
proj <- addImputeWeights(proj)
p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = 'GeneScoreMatrix',
  name = markerGenes,
  embedding = 'UMAP',
  imputeWeights = getImputeWeights(proj)
)
plotPDF(plotList = p,
        name = 'featureplot.pdf',
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)

## assign celltype annotation
proj@cellColData@listData[["cell_type"]] <- recode(
  proj@cellColData@listData[["Clusters"]],
  'C1' = 'Erythroid',
  'C2' = 'B cell',
  'C3' = 'Immature B cell',
  'C4' = 'Plasma cell',
  'C5' = 'T cell',
  'C6' = 'Thymocyte',
  'C7' = 'Endothelium',
  'C8' = 'Endothelium',
  'C9' = 'Macrophage',
  'C10' = 'Alveolar macrophage',
  'C11' = 'Macrophage',
  'C12' = 'Myeloid cell',
  'C13' = 'Proximal tubule',
  'C14' = 'Hepatocyte',
  'C15' = 'Foveolar cell',
  'C16' = 'Enterocyte',
  'C17' = 'Enterocyte',
  'C18' = 'Uterine epithelium',
  'C19' = 'Urothelial cell',
  'C20' = 'Basal-like',
  'C21' = 'Acinar cell',
  'C22' = 'Mucous neck cell',
  'C23' = 'Parietal cell',
  'C24' = 'Chief cell',
  'C25' = 'Pulmonary epithelium',
  'C26' = 'Pituitary stem cell',
  'C27' = 'Distal tubule and collecting duct',
  'C28' = 'Pancreatic secretory cell',
  'C29' = 'Cardiac muscle',
  'C30' = 'Skeletal muscle',
  'C31' = 'Smooth muscle cell',
  'C32' = 'Fibroblast',
  'C33' = 'Fibroblast',
  'C34' = 'Fibroblast',
  'C35' = 'Fibroblast',
  'C36' = 'Adipocyte',
  'C37' = 'Adrenal cortex',
  'C38' = 'Adrenal cortex',
  'C39' = 'Granulosa cell',
  'C40' = 'Granule cell',
  'C41' = 'Granule cell',
  'C42' = 'Neuron_Lin28b',
  'C43' = 'GABAergic neuron',
  'C44' = 'Endocrine cell',
  'C45' = 'Somatotrope',
  'C46' = 'Chromaffin cell',
  'C47' = 'Corticotrope/Gonadotrope',
  'C48' = 'Lactotrope',
  'C49' = 'Oligodendrocyte',
  'C50' = 'Astrocyte',
  'C51' = 'Oligodendrocyte',
  'C52' = 'Oligodendrocyte'
)

saveRDS(proj, file = file.path(ProjDir,'proj_Global.rds'))

## extract gene matrix
GeneScoreMat <- getMatrixFromProject(proj2,useMatrix = 'GeneScoreMatrix')
genemat <- GeneScoreMat@assays@data@listData$GeneScoreMatrix
rownames(genemat) <- as.character(GeneScoreMat@elementMetadata@listData[["name"]])



## ------------ sub-cluster for each lineage compartment ------------
## firstly, assign main cell type to one of 10 lineages: Endothelial, Epithelial, Erythroid, Hepatocyte, Immune, Muscle, Neuronal, Reproductive, Secretory, Stromal
## Then perform sub-clustering for each lineage and assign 'subanno' to every cell

## (the detailed code is omitted here for concision)
proj <- readRDS(file = file.path(ProjDir,'proj_Global_subanno.rds'))



## ------------ 5 add Peakset ------------
## Pseudo-bulk Replicates
proj2 <- addGroupCoverages(ArchRProj = proj, groupBy = 'subanno', minCells = 40, maxCells = 500)

## Calling peaks
proj2 <- addReproduciblePeakSet(
  ArchRProj = proj2,
  groupBy = "Clusters",
  pathToMacs2 = '~/anaconda3/bin/macs2',
  force = T
)
getPeakSet(proj2)
## Add Peak Matrix
proj2 <- addPeakMatrix(proj2)
getAvailableMatrices(proj2)

## Identify marker peaks
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj2,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
markerList


## plot Heatmap
options(bitmapType='cairo')
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  transpose = TRUE
)
pdf(file = file.path(picDir,'MarkersPeaksHeatmap.pdf'), width = 10, height = 7)
draw(heatmapPeaks,
     heatmap_legend_side = "bot", annotation_legend_side = "bot",
     # row_order = reorder_name,
     cluster_rows = T)
dev.off()



## ------------ 6 Motif Enrichment ------------
proj2 <- addMotifAnnotations(ArchRProj = proj2, motifSet = "cisbp", name = "Motif")
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj2,
  peakAnnotation = "Motif",
  cutOff = "Pval <= 0.01 & Log2FC >= 0.5"
)
enrichMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 12, height = 8, ArchRProj = projHeme5, addDOC = FALSE)



## ------------ 7 RNA Transfer ------------
## load RNA proj
seRNA <- readRDS(file = '/file/path/prefix/seRNA.rds')

proj2 <- addGeneIntegrationMatrix(
  ArchRProj = proj2,
  useMatrix = 'GeneScoreMatrix',
  matrixName = 'GeneIntegrationMatrix',
  reducedDims = 'IterativeLSI',
  seRNA = seRNA,
  addToArrow = TRUE,
  force = TRUE,
  # groupList = groupList,
  groupRNA = 'RNAlabel',
  nameCell = 'predictedCell',
  nameGroup = 'predictedGroup',
  nameScore = 'predictedScore'
)
saveRDS(proj2,file = file.path(ProjDir,'proj_Global_RNA.rds'))


