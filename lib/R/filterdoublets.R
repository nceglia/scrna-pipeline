library(Seurat)
library(SingleCellExperiment)

args = commandArgs(trailingOnly=TRUE)
seurat <- readRDS(args[1])
csv <- args[2]
sample_id <- args[3]
score <- args[4]

scores <- read.csv(csv,header=TRUE)

seurat <- AddMetaData(object=seurat, metadata=scores$Prediction, col.name='doublet')
seurat <- AddMetaData(object=seurat, metadata=scores$Score, col.name='doublet_score')
seurat <- AddMetaData(object=seurat, metadata=sample_id, col.name='sample_id')
seurat <- subset(seurat, subset = doublet_score < 1.0)
saveRDS(seurat,file=args[5])