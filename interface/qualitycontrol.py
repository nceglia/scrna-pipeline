from interface.tenxanalysis import TenxAnalysis

import subprocess
import os
import shutil
import uuid


class QualityControl(object):

    def __init__(self, tenx, sampleid, mouse=False):
        self.tenx = tenx
        self.sampleid = sampleid
        self.sce = os.path.join(tenx.path,"{0}.rdata".format(sampleid))
        self.qcdsce = os.path.join(tenx.path,"{0}_qcd.rdata".format(sampleid))
        self.cache = ".cache"
        if not os.path.exists(self.cache):
            os.makedirs(self.cache)
        self.construct = os.path.join(self.cache,"build_{}.R".format(sampleid))
        self.figures = os.path.join(self.cache,"figures_{}.R".format(sampleid))
        self.plots = os.path.join(self.tenx.path, "qc_figures")
        self.veryraw = os.path.join(self.tenx.path, "raw_build_{}.R".format(sampleid))

        output = open(self.construct,"w")
        output.write(script)
        output.close()

        output = open(self.figures,"w")
        output.write(figures)
        output.close()

        output = open(self.veryraw,"w")
        output.write(raw)
        output.close()
        if not os.path.exists(self.plots):
            os.makedirs(self.plots)

    def build(self, lsf=True):
        mat = self.tenx.filtered_matrices()
        cmd = ""
        if lsf:
            cmd = "/admin/lsf/10.1/linux3.10-glibc2.17-x86_64/bin/bsub -K -J cellassign -R rusage[mem=8] -n 4 -W 10 -o out -e err /opt/common/CentOS_7/singularity/singularity-3.0.1/bin/singularity exec --bind /admin --bind /opt --bind /work/shah/ /work/shah/images/base_scrna_r.img"
        subprocess.call(cmd.split() + ["Rscript", self.construct, mat, self.sce])

    def build_raw(self):
        mat = self.tenx.raw_matrices()
        print(" ".join(["Rscript", self.veryraw, mat, self.sce]))
        subprocess.call(["Rscript", self.veryraw, mat, self.sce])

    def plot(self):
        assert os.path.exists(self.sce), "SCE needs to be built before plotted."
        mat = self.tenx.filtered_matrices()
        print (" ".join(["Rscript", self.figures, self.sce, self.plots]))
        subprocess.call(["Rscript", self.figures, self.sce, self.plots])

    def run(self, mito=10):
        print("Running QC...")
        self.build()
        self.plot()

    def move(self, path):
        shutil.copyfile(self.sce, path)

    def sce(self):
        return SingleCellExperiment.fromRData(self.sce)

raw = """
library(scater)
library(SingleCellExperiment)
library(DropletUtils)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

sce <- read10xCounts(args[1])

saveRDS(sce, file=args[2])
print("Complete!")
"""


script = """
library(knitr)
library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(pheatmap)
library(ggrepel)
library(grid)
library(scran)
library(yaml)
library(scales)
library(feather)
library(limma)
library(Seurat)
library(scrna.utils)
library(scrna.sceutils)
library(cellassign)
library(cellassign.utils)
library(vdj.utils)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

args = commandArgs(trailingOnly=TRUE)
sce <- readRDS(args[1])
counts(sce) <- data.matrix(counts(sce))
logcounts(sce) <- data.matrix(logcounts(sce))
rownames(sce) <- rowData(sce)$Symbol
colnames(sce) <- colData(sce)$Barcode
sce$site <- args[2]
sce$patient <- args[3]
sce$sort <- args[4]
sce$therapy <- args[5]
sce$sample <- args[6]

mito_genes <- as.character(rowData(sce)$Symbol[str_detect(rowData(sce)$Symbol, "^MT\\-")]) %>% get_ensembl_id(sce)
ribo_genes <- as.character(rowData(sce)$Symbol[str_detect(rowData(sce)$Symbol, "^RP(L|S)")]) %>% unique %>% get_ensembl_id(sce)
sce <- calculateQCMetrics(sce, exprs_values = "counts", feature_controls =list(mitochondrial=mito_genes, ribosomal=ribo_genes))
sce <- filter_cells(sce, nmads = 3, type = "lower",log = TRUE, max_mito = 20, max_ribo = 60)

qclust <- quickCluster(sce, min.size = 100)
sce <- computeSumFactors(sce, clusters = qclust)
sce$size_factor <- sizeFactors(sce)
sce <- normalize(sce)

#Reduced dimensions
print("Running PCA")
sce <- runPCA(sce, ncomponents = 50, exprs_values = "logcounts")
print("Running TSNE")
sce <- runTSNE(sce, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)
print("Running UMAP")
sce <- runUMAP(sce, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)


saveRDS(sce, file=args[2])
print("Finished SCE Build")
"""


figures = """
library(scater)
library(SingleCellExperiment)
library(stringr)
library(scran)

args = commandArgs(trailingOnly=TRUE)

rdata <- readRDS(args[1])
sce <- as(rdata, 'SingleCellExperiment')

print("Generating Figures")

png(paste0(args[2],"/mito_v_features_by_counts.png"))
plotColData(sce, x = "total_features_by_counts", y = "pct_counts_mito")
dev.off(0)

png(paste0(args[2],"/total_counts_v_features_by_counts.png"))
plotColData(sce, x = "total_features_by_counts", y = "total_counts")
dev.off(0)

png(paste0(args[2],"/counts.png"))
hist(sce$total_counts, breaks=20, col='darkgoldenrod1',xlab='Total Counts')
dev.off(0)

png(paste0(args[2],"/umi.png"))
hist(sce$log10_total_counts, breaks=20, col='darkgoldenrod1',xlab='Log-total UMI count')
dev.off(0)

png(paste0(args[2],"/mito.png"))
hist(sce$pct_counts_mito, breaks=20, col='darkolivegreen4',xlab='Proportion of reads in mitochondrial genes')
dev.off(0)

png(paste0(args[2],"/ribo.png"))
hist(sce$pct_counts_ribo, breaks=20, col='firebrick4',xlab='Proportion of reads in ribosomal genes')

print("Finished Figures")
"""
