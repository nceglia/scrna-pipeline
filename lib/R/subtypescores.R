library(Seurat)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(scater)
library(SingleCellExperiment)

args = commandArgs(trailingOnly=TRUE)

seurat <- readRDS(args[1])
subtype <- args[2]
read_genes <- read.csv(args[3])
geneset_oi <- as.character(read_genes$genes)
subtype_figure <- args[4]

seurat <- AddModuleScore(object = seurat, features = geneset_oi, ctrl = 5, name = "subtype_score")
saveRDS(seurat, file="/work/shah/ceglian/tester.rdata")
subfig <- FeaturePlot(seurat,reduction="umapcelltype",features="subtype_score1") + ggtitle(paste0(subtype, " Score")) + theme(plot.title=element_text(face="bold")) + xlab("UMAP-1") + ylab("UMAP-2")

qdf <- data.frame(score=seurat$subtype_score1)
hist <- ggplot(qdf, aes(x=score)) + geom_histogram() + xlab('Score') + ylab("Cells") + ggtitle(paste0(subtype, " Score")) + theme(axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10), plot.title=element_text(face="bold"))


df <- data.frame(score=seurat$subtype_score1, sample=seurat$sample, batch=seurat$batch)
new_order <- with(df, reorder(sample, score, mean, na.rm=T))
df$sample <- factor(df$sample, levels = levels(new_order),ordered = TRUE)
top <- ggplot(as.data.frame(df), aes(x=batch, y=score, fill=sample)) + geom_boxplot()+  xlab("") + ylab(paste0(subtype, " Score"))  + ggtitle("Sample Summary") + theme(legend.position="none", plot.title=element_text(size=19, face = "bold"),axis.title.x=element_text(size=10, face = "bold"),axis.title.y=element_text(size=10, face = "bold"),axis.text.y = element_text(face="bold",size=10),axis.text.x = element_text(face="bold",size=10,angle = 90, hjust = 1),legend.title=element_text(size=0))

new_order <- with(df, reorder(batch, score, mean, na.rm=T))
df$batch <- factor(df$batch, levels = levels(new_order),ordered = TRUE)
bottom <- ggplot(as.data.frame(df), aes(x=sample, y=score, fill=batch)) + geom_boxplot() +  xlab("") + ylab(paste0(subtype, " Score"))  + ggtitle("Batch Summary") + theme(legend.position="none", plot.title=element_text(size=19, face = "bold"),axis.title.x=element_text(size=10, face = "bold"),axis.title.y=element_text(size=10, face = "bold"),axis.text.y = element_text(face="bold",size=10),axis.text.x = element_text(size=10,angle = 90, hjust = 1),legend.title=element_text(size=0))

xbottom <- cowplot::plot_grid(top, bottom, axis = 'tblr', ncol = 1, nrow=2, rel_widths = c(4, 10))

top <- plot_grid(subfig, hist, nrow = 1, ncol=2, align = "h") 
figure <- plot_grid(top, xbottom, nrow = 2, ncol=1, align = "v") 
ggsave(subtype_figure,figure,width=6,height=6, scale=2)
