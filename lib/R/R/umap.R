library(schex)
library(Seurat)
library(scater)

args = commandArgs(trailingOnly=TRUE)
annotated_seurat <- args[1]
umap   <- args[2]
column <- args[3]

seurat <- readRDS(annotated_seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, verbose = TRUE)
seurat <- RunUMAP(seurat, dims=1:50, verbose = TRUE)

Idents(seurat) <- "cell_type"
fig <- DimPlot(seurat,reduction="umap")
ggsave(umap, fig, width=2,height=2,scale=2)
