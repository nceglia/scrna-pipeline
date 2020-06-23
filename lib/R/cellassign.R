library(Seurat)
library(SingleCellExperiment)
library(scran)
library(reticulate)
library(cellassign)
library(tensorflow)
library(mltools)
library(data.table)
library(cowplot)
library(ggplot2)
library(tidyverse)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
merged_file      <- args[1]
marker_csv       <- args[2]
probabilities    <- args[3]
annotated_seurat <- args[4]
batch_figure     <- args[5]
batch_csv     <- args[6]

merged <- readRDS(merged_file)

DefaultAssay(object = merged) <- "RNA"
sce <- as.SingleCellExperiment(merged)
counts(sce) <- data.matrix(counts(sce))
sce <- computeSumFactors(sce)
sce$size_factor <- sizeFactors(sce)

rho <- read.csv(marker_csv)
rownames(rho) <- rho$X
rho <- rho[,-1]

sce_markers <- sce[rowSums(counts(sce)) > 1,]
genes <- intersect(rownames(sce_markers),rownames(rho))
sce_markers <- sce_markers[genes,]
sce_markers <- sce_markers[,colSums(counts(sce_markers))>1]
genes <- intersect(rownames(sce_markers),rownames(rho))
rho <- rho[genes,]
s <- sizeFactors(sce_markers)
rho <- data.matrix(rho)
# design <- data.frame(id=colnames(sce_markers),sample=sce_markers$sample)
# design <- as.data.table(design)
# design <- model.matrix(~sample, data = design)

#res <- cellassign(exprs_obj = sce_markers, X = design, s = s, marker_gene_info = rho, B = 20, shrinkage = TRUE, verbose = TRUE, rel_tol_em = 1e-4, num_runs=1)
res <- cellassign(exprs_obj = sce_markers, s = s, marker_gene_info = rho, B = 20, shrinkage = TRUE, verbose = TRUE, rel_tol_em = 1e-4, num_runs=1)
saveRDS(res, file=probabilities)

sce$cell_type <- "Other"
sce[,colnames(sce_markers)]$cell_type <- res$cell_type
merged$cell_type <- sce$cell_type

merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged, vars.to.regress = c("CC.Diff", "percent.mt", "nCount_RNA"))
merged <- RunPCA(merged, verbose = TRUE)
merged <- RunUMAP(merged, dims = 1:50)

merged <- FindNeighbors(object = merged)
merged <- FindClusters(object = merged, resolution=0.2)
markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
marker_table <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
marker_table <- as.data.frame(marker_table)
marker_table <- marker_table[c("p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","gene")]
write.csv(marker_table, file=batch_csv, row.names=FALSE)

saveRDS(merged, file=annotated_seurat)
Idents(merged) <- "cell_type"
df <- data.frame(cell_type=merged$cell_type, sample=merged$sample)
sample_summary <- summarise(group_by(df,cell_type,sample),count =n())
sample_summary <- data.frame(sample_summary)

fig <- ggplot(sample_summary, aes(factor(sample), count, fill = cell_type)) + geom_bar(stat="identity",position="fill") + xlab("Sample") + ylab("Cell Type") + theme(legend.title = element_blank(),plot.title=element_text(face="bold")) + ggtitle("Composition")

mplot <- DotPlot(merged, features = genes) + RotatedAxis() + ggtitle("Cell Type Markers")


Idents(merged) <- "cell_type"
ctumap <- DimPlot(merged,reduction="umap") + ggtitle("Cell Type") + xlab('UMAP-1') + ylab('UMAP-2')
Idents(merged) <- "sample"
sumap <- DimPlot(merged,reduction="umap")  + ggtitle("Sample") + xlab('UMAP-1') + ylab('UMAP-2')
Idents(merged) <- "seurat_clusters"
clust <- DimPlot(merged,reduction="umap")  + ggtitle("Cluster") + xlab('UMAP-1') + ylab('UMAP-2')


top <- plot_grid(ctumap, sumap, fig, nrow = 1, ncol=3, align = "h")
bottom <- plot_grid(clust, mplot, nrow = 1, ncol=2, rel_widths=c(1,2), align = "h")

figure <- plot_grid(top, bottom, nrow = 2, ncol=1, align = "v") + theme(text=element_text(size=12, family="Helvetica"))
ggsave(batch_figure,figure,width=8,height=4, scale=2)
