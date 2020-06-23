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
marker_csv <- args[5]


Idents(merged) <- "cell_type"
seurat <- subset(merged, idents = c(celltype))
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, verbose = TRUE)
seurat <- RunHarmony(seurat, group.by.vars = "batch", dims.use = 1:50, assay.use = "RNA")
seurat <- RunUMAP(seurat, dims = 1:50, reduction = "harmony", reduction.name = "umapcelltype", reduction.key = "umapcelltype_")
saveRDS(seurat, file=celltype_obj)

Idents(seurat) <- "sample"
sumap <- DimPlot(seurat,reduction="umapcelltype") + ggtitle("Sample") + xlab('UMAP-1') + ylab('UMAP-2')
Idents(seurat) <- "batch"
bumap <- DimPlot(seurat,reduction="umapcelltype") + ggtitle("Batch") + xlab('UMAP-1') + ylab('UMAP-2')

seurat <- FindNeighbors(object = seurat)
seurat <- FindClusters(object = seurat)
markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
marker_table <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
marker_table <- as.data.frame(marker_table)
marker_table <- marker_table[c("p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","gene")]
write.csv(marker_table, file=marker_csv, row.names=FALSE)

Idents(merged) <- "seurat.cluster"
cluster <- DimPlot(seurat,reduction="umapcelltype") + ggtitle("Cluster") + xlab('UMAP-1') + ylab('UMAP-2')

figure <- plot_grid(sumap, bumap, cluster, nrow = 1, align = "h") + theme(text=element_text(size=12, family="Helvetica"))
ggsave(batch_figure,figure,width=9,height=2.8, scale=1.5)

