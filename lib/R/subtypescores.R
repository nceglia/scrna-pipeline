library(Seurat)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(scater)
library(SingleCellExperiment)
library(schex)

args = commandArgs(trailingOnly=TRUE)

seurat <- readRDS(args[1])
subtype <- args[2]
read_genes <- read.csv(args[3])
geneset_oi <- c(as.character(read_genes$genes))
subtype_figure <- args[4]

seurat <- AddModuleScore(object = seurat, features = geneset_oi, name = "subtype_score")
binned <- make_hexbin(seurat, nbins=28, dimension_reduction = "umapcelltype")
subfig <- plot_hexbin_meta(binned, col="subtype_score1", action="median")+ ggtitle(paste0(subtype, " Score")) + theme(plot.title=element_text(face="bold")) + xlab("UMAP-1") + ylab("UMAP-2")

qdf <- data.frame(score=seurat$subtype_score1,batch=seurat$group)
hist <- ggplot(qdf, aes(x=score)) + geom_histogram() + xlab('Score') + ylab("Cells") + ggtitle(paste0(subtype, " Score")) + theme(axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10), plot.title=element_text(face="bold"))

Idents(seurat) <- "group"
binned <- make_hexbin(seurat, nbins=60, dimension_reduction = "umapcelltype")
dmp <- plot_hexbin_meta(binned, col="sample", action="majority")+ ggtitle("Group") + theme(plot.title=element_text(face="bold")) + xlab("UMAP-1") + ylab("UMAP-2") 

histbatch <- ggplot(qdf, aes(x=score)) + geom_histogram() + xlab('Score') + ylab("Cells") + ggtitle("") + theme(axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10), plot.title=element_text(face="bold")) + facet_wrap(~batch,ncol=1)

figure <- plot_grid(subfig, hist, dmp, histbatch, nrow = 2, ncol=2, align = "vh") 
ggsave(subtype_figure,figure,width=6,height=6, scale=2)
