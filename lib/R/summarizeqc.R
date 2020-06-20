library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(scater)
library(ggplot2)
library(cowplot)
library(ggpubr)

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
mito <- ggplot(qdf, aes(x=mito)) + geom_histogram() + xlab('Proportion of reads in mitochondrial genes') + ylab("Cells") + geom_vline(xintercept=25, linetype="dashed", color = "red") + ggtitle("Mitochondrial") + theme(axis.title.x = element_text(size = 8),axis.title.y = element_text(size = 8))
ribo <- ggplot(qdf, aes(x=ribo)) + geom_histogram() + xlab('Proportion of reads in ribosomal genes') + ylab("Cells") + ggtitle("Ribosomal") + theme(axis.title.x = element_text(size = 8),axis.title.y = element_text(size = 8))
doub <- ggplot(ddf, aes(x=doublet_score)) + geom_histogram() + xlab('Doublet Score') + ylab("Cells") + geom_vline(xintercept=0.25, linetype="dashed", color = "red") + ggtitle("Doublet Scores") + theme(axis.title.x = element_text(size = 8),axis.title.y = element_text(size = 8))

seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, dims=1:50)

# genes <- FeaturePlot(seurat,reduction="umap", features=markers$X)

sce_qc = as.SingleCellExperiment(seurat)
sce_qc = calculateQCMetrics(sce_qc, exprs_values = "counts")
umi = plotReducedDim(sce_qc,dimred="UMAP",colour_by="total_counts") + ggtitle("UMI Counts") + theme(legend.title = element_blank())
dumap = plotReducedDim(sce_qc,dimred="UMAP",colour_by="doublet_score") + ggtitle("Doublet Score") + theme(legend.title = element_blank())
ccumap = plotReducedDim(sce_qc,dimred="UMAP",colour_by="CC.Diff") + ggtitle("S - G2M Score") + theme(legend.title = element_blank())

figure <- plot_grid(mito, ribo, doub, umi, dumap, ccumap, nrow = 2, align = "vh")

ggsave(report,figure,width=9,height=6)
