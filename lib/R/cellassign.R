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
library(scater)
library(schex)

args = commandArgs(trailingOnly=TRUE)
merged_file      <- args[1]
marker_csv       <- args[2]
probabilities    <- args[3]
annotated_seurat <- args[4]
batch_figure     <- args[5]
batch_csv        <- args[6]

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

# if (length(unique(sce_markers$sample)) > 1) {
#     design <- data.frame(id=colnames(sce_markers),sample=sce_markers$sample)
#     design <- as.data.table(design)
#     design <- model.matrix(~sample, data = design)
#     res <- cellassign(exprs_obj = sce_markers, X = design, s = s, marker_gene_info = rho, B = 20, shrinkage = TRUE, verbose = TRUE, rel_tol_em = 1e-5, num_runs=3)
# } else {
#     res <- cellassign(exprs_obj = sce_markers, s = s, marker_gene_info = rho, B = 20, shrinkage = TRUE, verbose = TRUE, rel_tol_em = 1e-5, num_runs=3)
# }

# saveRDS(res, file=probabilities)

sce$cell_type <- "Cancer.cell"
# sce[,colnames(sce_markers)]$cell_type <- res$cell_type

# cancer <- sce[,sce$cell_type == "Cancer.cell"]
# cancer$cell_type <- as.character(lapply(cancer$group, function(x) {paste0(x,"-Cancer.Cell")}))
# sce[,colnames(cancer)]$cell_type <- cancer$cell_type

merged$cell_type <- sce$cell_type

merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged, vars.to.regress = c("CC.Diff", "percent.mt", "nCount_RNA"))
merged <- RunPCA(merged, verbose = TRUE)
merged <- RunUMAP(merged, dims = 1:50)

merged <- FindNeighbors(object = merged)
merged <- FindClusters(object = merged, resolution=0.4)

markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

marker_table <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
marker_table <- as.data.frame(marker_table)
marker_table <- marker_table[c("p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","gene")]   
marker_table <- subset(marker_table, !grepl("RPS", gene) )
marker_table <- subset(marker_table, !grepl("RPL", gene) )
write.csv(marker_table, file=batch_csv, row.names=FALSE)

saveRDS(merged, file=annotated_seurat)

df <- data.frame(cell_type=merged$cell_type, sample=merged$sample)
sample_summary <- summarise(group_by(df,cell_type,sample),count =n())
sample_summary <- data.frame(sample_summary)

fig <- ggplot(sample_summary, aes(factor(sample), count, fill = cell_type)) + geom_bar(stat="identity",position="fill") + xlab("Sample") + ylab("Cell Type") + theme(legend.title = element_blank(),plot.title=element_text(face="bold"),axis.text.x=element_text(size=9,angle=90)) + ggtitle("Composition") + scale_fill_manual(values=c("#9d65ff","#fa8419","#98e024","#f4005f","gray","#58d1eb","#98e024","#000000","navy"))

merged <- merged[,merged$cell_type!="Other"]
# binned <- make_hexbin(merged, nbins=70, dimension_reduction = "umap")
# ctumap <- plot_hexbin_meta(binned, col="cell_type", action="majority") + xlab("UMAP-1") + ylab("UMAP-2") + guides(fill = guide_legend(title = ""),override.aes = list(size=9)) + scale_fill_manual(values=c("#9d65ff","#fa8419","#98e024","#f4005f","gray","#58d1eb","#98e024","#000000","navy"))  + ggtitle("Cell Type") + theme(plot.title=element_text(size=15, face = "bold"),axis.title.x=element_text(size=9, face = "bold"),axis.title.y=element_text(size=9, face = "bold"), legend.text=element_text(size=9, face = "bold"),axis.text.y = element_text(face="bold",size=9),axis.text.x = element_text(face="bold",size=9))
Idents(merged) <- "cell_type"
ctumap <- DimPlot(merged, reduction="umap",label=TRUE)

binned <- make_hexbin(merged, nbins=70, dimension_reduction = "umap")
sumap <- plot_hexbin_meta(binned, col="sample", action="majority") + xlab("UMAP-1") + ylab("UMAP-2") + guides(fill = guide_legend(title = ""),override.aes = list(size=9))  + ggtitle("Sample") + theme(plot.title=element_text(size=15, face = "bold"),axis.title.x=element_text(size=9, face = "bold"),axis.title.y=element_text(size=9, face = "bold"), legend.text=element_text(size=9, face = "bold"),axis.text.y = element_text(face="bold",size=9),axis.text.x = element_text(face="bold",size=9))
# sumap <- DimPlot(merged, reduction="umap",label=TRUE)

clust <- plot_hexbin_meta(binned, col="seurat_clusters", action="majority") + xlab("UMAP-1") + ylab("UMAP-2") + guides(fill = guide_legend(title = ""),override.aes = list(size=9)) + ggtitle("Cluster") + theme(plot.title=element_text(size=15, face = "bold"),axis.title.x=element_text(size=9, face = "bold"),axis.title.y=element_text(size=9, face = "bold"), legend.text=element_text(size=9, face = "bold"),axis.text.y = element_text(face="bold",size=9),axis.text.x = element_text(face="bold",size=9))
# clust <- DimPlot(merged, reduction="umap",label=TRUE)

figure <- plot_grid(ctumap, fig, sumap, clust, nrow = 2, ncol=2, align = "vh")
ggsave(batch_figure,figure,width=7,height=7, scale=2)
