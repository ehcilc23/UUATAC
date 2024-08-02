library(ArchR)
library(ggplot2)
library(parallel)
library(dplyr)
addArchRThreads(threads = 16)
addArchRGenome("mm10")

########################################################################
FragmentsDir <- paste0("./")
inputFiles <- list.files(path = FragmentsDir, pattern = '*sample3k_fragmentuse.tsv.gz$', full.names = T)
inputFiles
sample <- c('UUATAC_sample3k')
names(inputFiles) <- sample
inputFiles

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 0,
  minFrags = 1,
  maxFrags = 1000000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  TileMatParams = list(tileSize = 500),
  threads = 16
)

############################## doublet
ArrowFiles <- c("UUATAC_sample3k.arrow")
doubScore <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1,
  threads=16
)

############################## create object
ArrowFiles <- c("UUATAC_sample3k.arrow")
proj1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "Output",
  copyArrows = FALSE
)
proj1
dim(proj1@cellColData)[1]

proj1 <- filterDoublets(proj1, filterRatio=1.7)
dim(proj1@cellColData)[1]
proj1@cellColData@listData[["nFrags"]] %>% median()
median(proj1@cellColData@listData[["TSSEnrichment"]])
median(proj1$DoubletEnrichment)
mean(proj1$DoubletEnrichment)

############################## Call peak
table(proj1$Sample)
proj1_use <- addGroupCoverages(ArchRProj = proj1, force = TRUE,
                               groupBy = "Sample",threads = 1,maxReplicates=8)
pathToMacs2 <- '/usr/local/bin/macs2'

proj1_use <- addReproduciblePeakSet(
  ArchRProj = proj1_use, force = TRUE,
  groupBy = "Sample", 
  pathToMacs2 = pathToMacs2)

getAvailableMatrices(proj1_use)
proj1_use <- addPeakMatrix(proj1_use)
save(proj1_use, file = 'proj_UUATAC_sample3k.rda')
