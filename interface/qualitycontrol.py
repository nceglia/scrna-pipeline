from interface.tenxanalysis import TenxAnalysis

from azure.common.client_factory import get_client_from_cli_profile
from azure.mgmt.compute import ComputeManagementClient
from azure.mgmt.network import NetworkManagementClient
from azure.storage.blob import BlockBlobService, PublicAccess

import subprocess
import os
import shutil
import uuid

aztok = open(".sas_token","r").read().strip()

class QualityControl(object):

    def __init__(self, tenx, sampleid, mouse=False):
        self.tenx = tenx
        self.sampleid = sampleid
        self.sce = os.path.join(tenx.path,"{0}.rdata".format(sampleid))
        self.qcdsce = os.path.join(tenx.path,"{0}_qcd.rdata".format(sampleid))
        self.cache = ".cache"
        if not os.path.exists(self.cache):
            os.makedirs(self.cache)
        self.script = os.path.join(self.cache,"qc_{}.R".format(sampleid))
        self.construct = os.path.join(self.cache,"build_{}.R".format(sampleid))
        self.figures = os.path.join(self.cache,"figures_{}.R".format(sampleid))
        self.plots = os.path.join(self.tenx.path, "qc_figures")
        self.veryraw = os.path.join(self.tenx.path, "raw_build_{}.R".format(sampleid))

        output = open(self.script,"w")
        output.write(filter)
        output.close()

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
        self.storage_account = "scrnadata"
        try:
            version = self.tenx.detected_version
        except Exception as e:
            version = "v3"
        if mouse:
            self.container = "rdatamouse{}".format(version)
            self.rawcontainer = "rdatarawmouse{}".format(version)
        else:
            self.container = "rdata{}".format(version)
            self.rawcontainer = "rdataraw{}".format(version)
        self.block_blob_service = BlockBlobService(account_name='scrnadata', sas_token=aztok)

    def filter(self, mito=10):
        assert os.path.exists(self.sce), "SCE needs to be built before filtered."
        print (" ".join(["Rscript", self.script, self.sce, self.qcdsce, str(mito)]))
        subprocess.call(["Rscript", self.script, self.sce, self.qcdsce, str(mito)])

    def build(self):
        mat = self.tenx.filtered_matrices()
        print(" ".join(["Rscript", self.construct, mat, self.sce]))
        subprocess.call(["Rscript", self.construct, mat, self.sce])

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
        self.filter(mito=mito)
        self.plot()

    def move(self, path):
        shutil.copyfile(self.sce, path)

    def sce(self):
        return SingleCellExperiment.fromRData(self.sce)

    def upload(self):
        print ("Uploading QCD ", self.qcdsce, self.container)
        self.block_blob_service.create_blob_from_path(self.container, "{}.rdata".format(self.sampleid), self.qcdsce)

    def upload_raw(self):
        print ("Uploading RAW ", self.sce, self.container)
        self.block_blob_service.create_blob_from_path(self.rawcontainer, "{}.rdata".format(self.sampleid), self.sce)



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
library(scater)
library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(DropletUtils)
library(stringr)
library(scran)
library(annotables)

args = commandArgs(trailingOnly=TRUE)

sce <- read10xCounts(args[1])

rowData(sce)$ensembl_gene_id <- rownames(sce)

at <- annotables::grch38 # If GRCh38 (I think this is default?)
at <- at[!duplicated(at$ensgene),]
rd <- as.data.frame(rowData(sce))
rd$ensgene <- rd$ID
rd <- dplyr::left_join(rd, at, by = "ensgene")
rowData(sce) <- rd

print("Calculating Size Factors")
# Calculate size factors
sce_result <- tryCatch({scran::computeSumFactors(sce)},error=function(e) {NULL})
poolsize <- 100;
while (is.null(sce_result) && poolsize >= 0) {
    sce_result <- tryCatch({scran::computeSumFactors(sce,sizes=poolsize)},error=function(e) {NULL})
    if (is.null(sce_result)) {
      poolsize <- poolsize - 10;
    }
}

# Compute log normal expression values
sce <- normalize(sce)


# Get Mitochondrial genes for QC:
mt_genes <- as.character(rowData(sce)$Symbol[str_detect(rowData(sce)$Symbol, "^MT\\\-")])
ribo_genes <- as.character(rowData(sce)$Symbol[str_detect(rowData(sce)$Symbol, "^RP(L|S)")])
# Calculate QC metrics
print("Calculating QC Metrics")

rownames(sce) <- rowData(sce)$Symbol

sce <- calculateQCMetrics(sce, exprs_values = "counts", feature_controls =
                          list(mito=mt_genes, ribo=ribo_genes))
#Reduced dimensions
print("Running PCA")
sce <- runPCA(sce, ntop = 1000, ncomponents = 50, exprs_values = "logcounts")
print("Running TSNE")
sce <- runTSNE(sce, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)
print("Running UMAP")
sce <- runUMAP(sce, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)


saveRDS(sce, file=args[2])
print("Finished SCE Build")
"""


filter = """
library(scater)
library(SingleCellExperiment)
library(stringr)
library(scran)

args = commandArgs(trailingOnly=TRUE)

rdata <- readRDS(args[1])
sce <- as(rdata, 'SingleCellExperiment')

print("Filtering")
cells_to_keep <- sce$pct_counts_mito < args[3]
table_cells_to_keep <- table(cells_to_keep)
sce <- sce[,cells_to_keep]
# summ <- summary(sce$total_counts)
# thresh <- summ[[2]]
# keep <- sce$total_counts > thresh
# sce <- sce[,keep]

saveRDS(sce, file=args[2])
print("Finished Filtering")
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
