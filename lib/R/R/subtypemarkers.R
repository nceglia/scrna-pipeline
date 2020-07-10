library(Seurat)
library(scater)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
seurat <- readRDS(args[1])
marker_csv <- args[2]

seurat <- FindNeighbors(object = seurat)
seurat <- FindClusters(object = seurat)
markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
marker_table <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
marker_table <- as.data.frame(marker_table)
write.csv(marker_table, file = marker_csv)
