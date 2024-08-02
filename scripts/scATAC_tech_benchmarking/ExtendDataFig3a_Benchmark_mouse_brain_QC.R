library(ArchR)
library(ggplot2)
library(parallel)
library(dplyr)
addArchRThreads(threads = 16)
addArchRGenome("mm10")


########################################################################
FragmentsDir <- paste0("./")
inputFiles <- list.files(path = FragmentsDir, pattern = '*sort.bed.gz$', full.names = T)
inputFiles
sample <- c('UUATAC')
names(inputFiles) <- sample
inputFiles

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 0,
  minFrags = 10,
  maxFrags = 1000000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  TileMatParams = list(tileSize = 500), 
  threads = 16 
)

ArrowFiles <- c("UUATAC.arrow")
############################## Create an ArchRproject
proj1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "Output",
  copyArrows = FALSE,
)
proj1
dim(proj1@cellColData)[1]
save(proj1, file = 'proj_raw_UUATAC.rds')
df <- getCellColData(ArchRProj = proj_tmp, select = c('nFrags','TSSEnrichment')) %>% as.data.frame()
write.table(df,file = 'metadata_UUATAC.txt',sep = '\t',quote = F)


############################## QC plot
tss_cutoff = 6; uq_cutoff = 3000
idxPass <- which( (proj1$TSSEnrichment >= tss_cutoff) & (proj1$nFrags >= uq_cutoff) )
length(idxPass)

cellsPass <- proj1$cellNames[idxPass]
proj1 = proj1[cellsPass, ]
dim(proj1@cellColData)[1]

uq_value = median(proj1$nFrags); tss_value = median(proj1$TSSEnrichment); cell_value = length(idxPass)
pdf(file = 'UUATAC_qc.pdf', width = 5, height = 5, onefile = F)
ggPoint(
  x = df[,1], y = df[,2], colorDensity = TRUE, rastr = T,
  size = 1, continuousSet = "sambaNight",
  title = paste0("UUATAC Mouse Brain\nNumber of nuclei passing QC: ", cell_value, "\nMedian UQ: ", uq_value, "; Median TSS: ", tss_value, "\nCut-off: TSS > ", tss_cutoff, "; UQ > ", uq_cutoff),
  xlabel = "Log10 (Unique Fragments)", ylabel = "TSS Enrichment", 
) +
  geom_hline(yintercept = tss_cutoff, lty = "dashed") + 
  geom_vline(xintercept = log10(uq_cutoff), lty = "dashed") +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=10), plot.title=element_text(size=10), legend.position = 'right', legend.direction = "vertical")
dev.off()

cellsPass_name <- as.data.frame(cellsPass)
write.csv(cellsPass_name, file = 'cellspass_QC_UUATAC.csv', row.names = F)
