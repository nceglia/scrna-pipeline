library(Seurat)
library(harmony)
library(scater)
library(tidyverse)
library(cowplot)
library(ggplot2)

options(future.globals.maxSize= 2097152000)

args = commandArgs(trailingOnly=TRUE)

merged <- readRDS(args[1])
integrated_object_file <- args[2]
project_figure <- args[3]
celltype_csv <- args[4]

merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged, verbose = TRUE)

integrated <- RunHarmony(merged, group.by.vars="batch", dims.use = 1:50)
integrated <- RunPCA(integrated, verbose = TRUE)
integrated <- RunUMAP(integrated, dims = 1:50, reduction = "harmony", reduction.name = "umapharmony", reduction.key = "umapharmony_")
saveRDS(integrated, file=integrated_object_file)

Idents(integrated) <- "cell_type"
markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
marker_table <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
marker_table <- as.data.frame(marker_table)
marker_table <- marker_table[c("p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","gene")]
write.csv(marker_table, file=celltype_csv, row.names=FALSE)


Idents(integrated) <- "cell_type"
int <- DimPlot(integrated,reduction="umapharmony") + ggtitle("Integrated Cell Type") + xlab('UMAP-1') + ylab('UMAP-2')

Idents(integrated) <- "sample"
sam <- DimPlot(integrated,reduction="umapharmony") + ggtitle("Integrated Sample") + xlab('UMAP-1') + ylab('UMAP-2')

df <- data.frame(cell_type=merged$cell_type, batch=merged$batch)
batch_summary <- summarise(group_by(df,cell_type,batch),count =n())
batch_summary <- data.frame(batch_summary)

ct <- ggplot(batch_summary, aes(factor(batch), count, fill = cell_type)) + geom_bar(stat="identity",position="fill") + xlab("Batch") + ylab("Cell Type") + theme(legend.title = element_blank(),axis.text.x=element_text(angle=90)) + ggtitle("Celltypes")

totalct <- ggplot(df, aes(cell_type))+ geom_bar() + ggtitle("Cell Types") + xlab('') + ylab('Cells')+theme(legend.title = element_blank(),axis.text.x=element_text(angle=90))

figure <- plot_grid(int, sam, ct, totalct, nrow = 2, ncol=2, align = "h") + theme(text=element_text(size=12, family="Helvetica"))
ggsave(project_figure,figure,width=6,height=4, scale=2)

