library(Seurat)
library(scater)
library(tidyverse)
library(harmony)
library(cowplot)
library(ggplot2)
library(schex)
library(SingleCellExperiment)

args = commandArgs(trailingOnly=TRUE)
merged <- readRDS(args[1])
celltype <- args[2]
celltype_obj <- args[3]
batch_figure <- args[4]
marker_csv <- args[5]
sce_obj <- args[6]

Idents(merged) <- "cell_type"
seurat <- merged[,merged$cell_type==celltype]

seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, verbose = TRUE)
seurat <- RunHarmony(seurat, group.by.vars = "batch", dims.use = 1:50, assay.use = "RNA")
seurat <- RunUMAP(seurat, dims = 1:50, reduction = "harmony", reduction.name = "umapcelltype", reduction.key = "umapcelltype_")
saveRDS(seurat, file=celltype_obj)

binned <- make_hexbin(seurat, nbins=35, dimension_reduction = "umapcelltype")
sumap <- plot_hexbin_meta(binned, col="seurat_clusters", action="majority") + theme(legend.title=element_text(size=0,face = "bold"),legend.text=element_text(size=7, face = "bold"),plot.title=element_text(size=12, face = "bold")) + ggtitle("Cluster") + xlab('UMAP-1') + ylab('UMAP-2')

bumap <-  plot_hexbin_meta(binned, col="batch", action="majority") + theme(legend.title=element_text(size=0,face = "bold"),legend.text=element_text(size=7, face = "bold"),plot.title=element_text(size=12, face = "bold")) + ggtitle("Batch") + xlab('UMAP-1') + ylab('UMAP-2')

seurat <- FindNeighbors(object = seurat)
seurat <- FindClusters(object = seurat, resolution=0.1)
sce <- as.SingleCellExperiment(seurat)
saveRDS(sce, file=sce_obj)

markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
marker_table <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
marker_table <- as.data.frame(marker_table)
marker_table <- marker_table[c("p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","gene")]
marker_table <- subset(marker_table, !grepl("RPS", gene) )
marker_table <- subset(marker_table, !grepl("RPL", gene) )
write.csv(marker_table, file=marker_csv, row.names=FALSE)

binned <- make_hexbin(seurat, nbins=60, dimension_reduction = "umapcelltype")
cluster <-  plot_hexbin_meta(binned, col="sample", action="majority") + theme(legend.title=element_text(size=0,face = "bold"),legend.text=element_text(size=7, face = "bold"),plot.title=element_text(size=12, face = "bold")) + ggtitle("Sample") + xlab('UMAP-1') + ylab('UMAP-2')
df <- data.frame(sample=seurat$seurat_clusters)

totalsamp <- ggplot(df, aes(sample, fill=sample))+ geom_bar() + ggtitle("Cluster") + xlab('') + ylab('Cells')+theme(legend.position="none",axis.text.x=element_text(angle=90,size=7))

figure <- plot_grid(sumap, bumap, cluster, totalsamp, nrow = 2, ncol=2, align = "vh")
ggsave(batch_figure,figure,width=6,height=6, scale=2.0)

