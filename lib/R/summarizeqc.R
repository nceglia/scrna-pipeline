library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(scater)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(grid)
library(schex)

args = commandArgs(trailingOnly=TRUE)
seurat   <- readRDS(args[1])
sce      <- readRDS(args[2])
doublets <- read.csv(args[3],stringsAsFactors=FALSE)
markers  <- read.csv(args[4],stringsAsFactors=FALSE)
sample   <- args[5]
mthresh  <- args[6]
dthresh  <- args[7]
report   <- args[8]


saveRDS(markers, file="/work/shah/ceglian/mm.rdata")
mitochondrial <- as.character(rownames(sce)[str_detect(rownames(sce), "^MT\\-")])
ribosomal <- as.character(rownames(sce)[str_detect(rownames(sce), "^RP(L|S)")])

sce <- calculateQCMetrics(sce, exprs_values = "counts", feature_controls = list(mitochondrial=mitochondrial, ribosomal=ribosomal))

qdf <- data.frame(mito=sce$pct_counts_mitochondrial, ribo=sce$pct_counts_ribosomal)
ddf <- data.frame(doublet_score=doublets$Score)
saveRDS(df, file="/work/shah/ceglian/df.rdata")
mito <- ggplot(qdf, aes(x=mito)) + geom_histogram() + xlab('Proportion of reads in mitochondrial genes') + ylab("Cells") + geom_vline(xintercept=25, linetype="dashed", color = "red") + ggtitle("Mitochondrial") + theme(axis.title.x = element_text(size = 8),axis.title.y = element_text(size = 8), plot.title=element_text(face="bold"))
ribo <- ggplot(qdf, aes(x=ribo)) + geom_histogram() + xlab('Proportion of reads in ribosomal genes') + ylab("Cells") + ggtitle("Ribosomal") + theme(axis.title.x = element_text(size = 8),axis.title.y = element_text(size = 8), plot.title=element_text(face="bold"))
doub <- ggplot(ddf, aes(x=doublet_score)) + geom_histogram() + xlab('Doublet Score') + ylab("Cells") + geom_vline(xintercept=0.25, linetype="dashed", color = "red") + ggtitle("Doublet Scores") + theme(axis.title.x = element_text(size = 8),axis.title.y = element_text(size = 8), plot.title=element_text(face="bold"))

seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, dims=1:50)

sce_qc = as.SingleCellExperiment(seurat)
sce_qc$pct_counts_mitochondrial <- sce[,colnames(sce_qc)]$pct_counts_mitochondrial

sce_qc = calculateQCMetrics(sce_qc, exprs_values = "counts")

binned <- make_hexbin(sce_qc, nbins=30, dimension_reduction = "UMAP")

umi = plot_hexbin_meta(binned, col="log10_total_counts", action="median")  + ggtitle("UMI Log Counts") + theme(legend.title = element_blank())
features = plot_hexbin_meta(binned, col="log10_total_features_by_counts", action="median") + ggtitle("Total Features by Counts") + theme(legend.title = element_blank())
dumap = plot_hexbin_meta(binned, col="doublet_score", action="median") + ggtitle("Doublet Score") + theme(legend.title = element_blank())
ccumap = plot_hexbin_meta(binned, col="Phase", action="majority") + ggtitle("Cell Cycle Phase") + theme(legend.title = element_blank()) +scale_fill_manual(values=c("#9d65ff","#fa8419","#98e024"))
mitomap = plot_hexbin_meta(binned, col="pct_counts_mitochondrial", action="median") + ggtitle("Mito Percentage") + theme(legend.title = element_blank())

summary_stats = data.frame(Metric=c("Cells","Median Counts per Cell","Median Genes per Cell"), Value=c(length(colnames(seurat)),median(seurat$nCount_RNA),median(seurat$nFeature_RNA)))
stat <- ggtexttable(summary_stats, rows = NULL, theme = ttheme(base_style="mOrange", base_size=20))

cvm <- plotColData(sce, x = "log10_total_features_by_counts", y="pct_counts_mitochondrial")  + ggtitle("Feature vs. Percent Mito")


figure <- plot_grid(stat, mito, ribo, doub, umi, dumap, ccumap, mitomap, cvm, features, nrow = 5, ncol=2, align = "vh")

ggsave(report, figure, width=8,height=16, scale=1.5)
