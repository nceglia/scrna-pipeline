library(Seurat)
library(scater)
library(tidyverse)
library(harmony)
library(cowplot)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
merged <- readRDS(args[1])
celltype <- args[2]
celltype_obj <- args[3]
batch_figure <- args[4]
Idents(merged) <- "cell_type"
seurat <- subset(merged, idents = celltype)
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, verbose = TRUE)
seurat <- RunHarmony(seurat, group.by.vars = "batch", dims.use = 1:50, assay.use = "RNA")
seurat <- RunUMAP(seurat, dims = 1:50, reduction = "harmony", reduction.name = "umapharmony", reduction.key = "umapharmony_")
saveRDS(seurat, file=celltype_obj)

Idents(merged) <- "sample"
sumap <- DimPlot(merged,reduction="umapharmony")
Idents(merged) <- "batch"
bumap <- DimPlot(merged,reduction="umapharmony")

figure <- plot_grid(sumap, bumap, nrow = 1, align = "h")
ggsave(batch_figure,figure,width=8,height=4)
