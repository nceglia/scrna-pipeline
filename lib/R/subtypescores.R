library(Seurat)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(scater)
library(SingleCellExperiment)
library(schex)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)

seurat <- readRDS(args[1])
subtype <- args[2]
read_genes <- read.csv(args[3])
geneset_oi <- c(as.character(read_genes$genes))
subtype_figure <- args[4]

seurat <- AddModuleScore(object = seurat, features = geneset_oi, name = "subtype_score")
binned <- make_hexbin(seurat, nbins=28, dimension_reduction = "umapcelltype")
subfig <- plot_hexbin_meta(binned, col="subtype_score1", action="median")+ ggtitle(paste0(subtype, " Score")) + theme(plot.title=element_text(face="bold")) + xlab("UMAP-1") + ylab("UMAP-2")

qdf <- data.frame(score=seurat$subtype_score1,group=seurat$group,batch=seurat$batch)
hist <- ggplot(qdf, aes(x=score)) + geom_histogram() + xlab('Score') + ylab("Cells") + ggtitle(paste0(subtype, " Score")) + theme(axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10), plot.title=element_text(face="bold"))

sce <- as.SingleCellExperiment(seurat)
sce$subtype <- "True"
sce[,sce$subtype_score1 < 0.5]$subtype <- "False"
seurat$subtype <- sce$subtype

df <- data.frame(cell_type=seurat$subtype, group=seurat$group)
group_summary <- summarise(group_by(df,cell_type,group),count =n())
group_summary <- data.frame(group_summary)
colormap <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(seurat$group)))
binned <- make_hexbin(seurat, nbins=60, dimension_reduction = "umapcelltype")
dmp <- plot_hexbin_meta(binned, col="group", action="majority")+ ggtitle("Group") + theme(plot.title=element_text(face="bold")) + xlab("UMAP-1") + ylab("UMAP-2") + scale_fill_manual(values=colormap)

histbatch <- ggplot(qdf, aes(x=score)) + geom_histogram() + xlab('Score') + ylab("Cells") + ggtitle("") + theme(axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10), plot.title=element_text(face="bold")) + facet_wrap(~group,ncol=1) + xlim(0,3)

colormap <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(seurat$group)))
hist <- ggplot(group_summary, aes(factor(group), count, fill = cell_type)) + geom_bar(stat="identity",position="fill") + xlab("Group") + ylab("") + theme(axis.text.x=element_text(angle=90)) + ggtitle("Subtype by Group")+ scale_fill_manual(values=colormap)

figure <- plot_grid(subfig, histbatch, dmp, hist, nrow = 2, ncol=2, align = "vh") 
ggsave(subtype_figure,figure,width=6,height=6, scale=2)
