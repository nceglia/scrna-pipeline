library(Seurat)
library(tidyverse)
library(scater)

args = commandArgs(trailingOnly=TRUE)
seurat <- readRDS(args[1])
object_file <- args[2]

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seurat <- CellCycleScoring(seurat, s.genes, g2m.genes)
seurat$CC.Diff <- seurat$S.Score - seurat$G2M.Score
saveRDS(seurat, object_file)
