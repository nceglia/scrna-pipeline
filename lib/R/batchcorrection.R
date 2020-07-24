library(Seurat)
library(harmony)
library(scater)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(schex)

options(future.globals.maxSize= 2097152000)

args = commandArgs(trailingOnly=TRUE)

merged <- readRDS(args[1])
integrated_object_file <- args[2]
project_figure <- args[3]
celltype_csv <- args[4]
merged <- merged[,merged$cell_type!="Other"]
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged, verbose = TRUE)
merged <- RunUMAP(merged, dims = 1:50)

integrated <- RunHarmony(merged, group.by.vars="batch", dims.use = 1:50)
integrated <- RunPCA(integrated, verbose = TRUE)
integrated <- RunUMAP(integrated, dims = 1:50, reduction = "harmony", reduction.name = "umapharmony", reduction.key = "umapharmony_")
saveRDS(integrated, file=integrated_object_file)

integrated <- FindNeighbors(object = integrated)
integrated <- FindClusters(object = integrated, resolution=0.4)
Idents(integrated) <- "group"
markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
marker_table <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
marker_table <- as.data.frame(marker_table)
marker_table <- marker_table[c("p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","gene")]
marker_table <- subset(marker_table, !grepl("RPS", gene) )
marker_table <- subset(marker_table, !grepl("RPL", gene) )
write.csv(marker_table, file=celltype_csv, row.names=FALSE)


integrated <- integrated[,integrated$cell_type!="Other"]

binned_int <- make_hexbin(integrated, nbins=40, dimension_reduction = "umapharmony")
binned_m <- make_hexbin(integrated, nbins=40, dimension_reduction = "umap")
# binned_merged <- make_hexbin(merged, nbins=35, dimension_reduction = "umap")
mint <- plot_hexbin_meta(binned_m, col="sample", action="majority") + xlab("UMAP-1") + ylab("UMAP-2") + guides(fill = guide_legend(title = ""),override.aes = list(size=9)) + ggtitle("Sample") + theme(plot.title=element_text(size=14, face = "bold"),axis.title.x=element_text(size=9, face = "bold"),axis.title.y=element_text(size=9, face = "bold"), legend.text=element_text(size=9, face = "bold"),axis.text.y = element_text(face="bold",size=9),axis.text.x = element_text(face="bold",size=9),legend.title=element_text(size=9))


msam <- plot_hexbin_meta(binned_int, col="sample", action="majority") + xlab("UMAP-1") + ylab("UMAP-2") + guides(fill = guide_legend(title = ""),override.aes = list(size=4))  + ggtitle("Sample") + theme(plot.title=element_text(size=14, face = "bold"),axis.title.x=element_text(size=9, face = "bold"),axis.title.y=element_text(size=4, face = "bold"), legend.text=element_text(size=4, face = "bold"),axis.text.y = element_text(face="bold",size=4),axis.text.x = element_text(face="bold",size=4),legend.title=element_text(size=4))


# int <- plot_hexbin_meta(binned_int, col="cell_type", action="majority") + xlab("UMAP-1") + ylab("UMAP-2") + guides(fill = guide_legend(title = ""),override.aes = list(size=9)) + ggtitle("Cell Type") + theme(plot.title=element_text(size=14, face = "bold"),axis.title.x=element_text(size=9, face = "bold"),axis.title.y=element_text(size=9, face = "bold"), legend.text=element_text(size=9, face = "bold"),axis.text.y = element_text(face="bold",size=9),axis.text.x = element_text(face="bold",size=9),legend.title=element_text(size=9))

Idents(integrated) <- "sample"
int <- DimPlot(integrated, reduction="umapharmony",label=FALSE)
Idents(integrated) <- "sample"
mergedint <- DimPlot(integrated, reduction="umap",label=FALSE)

sam <- plot_hexbin_meta(binned_int, col="batch", action="majority") + xlab("UMAP-1") + ylab("UMAP-2") + guides(fill = guide_legend(title = ""),override.aes = list(size=4)) + ggtitle("Batch Integrated") + theme(plot.title=element_text(size=14, face = "bold"),axis.title.x=element_text(size=9, face = "bold"),axis.title.y=element_text(size=4, face = "bold"), legend.text=element_text(size=4, face = "bold"),axis.text.y = element_text(face="bold",size=4),axis.text.x = element_text(face="bold",size=4),legend.title=element_text(size=4))

mint <- plot_hexbin_meta(binned_m, col="batch", action="majority") + xlab("UMAP-1") + ylab("UMAP-2") + guides(fill = guide_legend(title = ""),override.aes = list(size=4)) + ggtitle("Batch Merged") + theme(plot.title=element_text(size=14, face = "bold"),axis.title.x=element_text(size=9, face = "bold"),axis.title.y=element_text(size=4, face = "bold"), legend.text=element_text(size=4, face = "bold"),axis.text.y = element_text(face="bold",size=4),axis.text.x = element_text(face="bold",size=4),legend.title=element_text(size=4))

df <- data.frame(cell_type=merged$cell_type, batch=merged$batch)
batch_summary <- summarise(group_by(df,cell_type,batch),count =n())
batch_summary <- data.frame(batch_summary)

ct <- ggplot(batch_summary, aes(factor(batch), count, fill = cell_type)) + geom_bar(stat="identity",position="fill") + xlab("Batch") + ylab("Cell Type") + theme(legend.position = "none",axis.text.x=element_text(angle=90)) + ggtitle("Celltypes")+ scale_fill_manual(values=c("#9d65ff","#fa8419","#98e024","#f4005f","gray","#58d1eb","#98e024","#000000","navy"))

totalct <- ggplot(df, aes(cell_type, fill=cell_type))+ geom_bar() + ggtitle("Cell Types") + xlab('') + ylab('Cells')+theme(legend.position="none",axis.text.x=element_text(angle=90))+ scale_fill_manual(values=c("#9d65ff","#fa8419","#98e024","#f4005f","gray","#58d1eb","#98e024","#000000","navy"))

figure <- plot_grid(int, mergedint, sam, mint, ct, totalct, nrow = 3, ncol=2, align = "vh")
ggsave(project_figure,figure,width=6,height=8, scale=2)
