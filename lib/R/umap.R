library(schex)
library(Seurat)

args = commandArgs(trailingOnly=TRUE)
annotated_seurat <- args[1]
umap   <- args[2]

seurat <- readRDS(annotated_seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, verbose = TRUE)
seurat <- RunUMAP(seurat, dims=1:50, verbose = TRUE)

binned <- make_hexbin(seurat, nbins=100, dimension_reduction = "UMAP")

fig <- plot_hexbin_meta(binned, col="cell_type", action="majority") + guides(fill = guide_legend(title = "")) + xlab("UMAP-1") + ylab("UMAP-2")

ggsave(umap, fig, scale=0.2, width=3,height=3)
