library(SingleCellExperiment)
library(Seurat)

args = commandArgs(trailingOnly=TRUE)
paths <- read.csv(args[1], sep="\t", header=TRUE)
object_file <- args[2]

i <- 0
rdata <- as.character(paths$path)
samples <- as.character(paths$sample)
seurats <- vector("list", length(rdata)-1)

init_seurat <- readRDS(rdata[[1]])
init_seurat <- RenameCells(object = init_seurat, add.cell.id = samples[1])
DefaultAssay(object = init_seurat) <- "RNA"
init_seurat$sample <-samples[1]
for (i in 2:length(rdata)) {
  seurat <- readRDS(rdata[i])
  seurat <- RenameCells(object = seurat, add.cell.id = samples[i])
  DefaultAssay(object = seurat) <- "RNA"
  seurat$sample <- samples[i]
  seurats[[i-1]] <- seurat
}
merged <- merge(init_seurat, y = seurats, project = "RNASCP")
saveRDS(merged, file=object_file)
