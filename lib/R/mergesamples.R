library(SingleCellExperiment)
library(Seurat)

args = commandArgs(trailingOnly=TRUE)

object_file <- args[2]

paths <- read.csv(args[1], sep="\t", header=TRUE)

i <- 0
rdata <- as.character(paths$path)
samples <- as.character(paths$sample)
groups <- as.character(paths$group)
seurats <- vector("list", length(rdata)-1)
message("loading init.")
init_seurat <- readRDS(rdata[[1]])
init_seurat <- RenameCells(object = init_seurat, add.cell.id = samples[1])
DefaultAssay(object = init_seurat) <- "RNA"
init_seurat$sample <- samples[1]
init_seurat$group <- groups[[1]]
message("init loaded.")
if (length(samples) == 1) {
  merged <- init_seurat
} else {
  for (i in 2:length(rdata)) {
    message(samples[i])
    seurat <- readRDS(rdata[i])
    seurat <- RenameCells(object = seurat, add.cell.id = samples[i])
    DefaultAssay(object = seurat) <- "RNA"
    seurat$sample <- samples[i]
    seurat$group <- groups[i]
    seurats[[i-1]] <- seurat
  }
  merged <- merge(init_seurat, y = seurats, project = "RNASCP")
}
message("finished samples.")
if (file.exists(args[3])) {
  phase <- merged$Phase
  annotated <- readRDS(args[3])
  common <- intersect(colnames(merged),colnames(annotated))
  merged <- annotated[,common]
  #merged$Phase <- "NA"  
}
saveRDS(merged, file=object_file)
