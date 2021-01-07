library(Seurat)
library(SingleCellExperiment)
library(scran)
library(reticulate)
library(cellassign)
library(tensorflow)
library(mltools)
library(data.table)
library(tidyverse)
library(dplyr)
library(scater)
library(scater)

args = commandArgs(trailingOnly=TRUE)
batch_merged_h5ad <- args[1]
marker_csv        <- args[2]
celltype_csv      <- args[3]

merged <- ReadH5AD(batch_merged_h5ad)
sce <- as.SingleCellExperiment(merged)
if ("cell_type" %in% colnames(colData(sce)))  {
    print("Loaded annotated seurat.")
} else {
    DefaultAssay(object = merged) <- "RNA"

    sce <- as.SingleCellExperiment(merged)
    symbols <- lapply(rownames(sce), function(x) { toupper(x) })
    rownames(sce) <- symbols
    rowData(sce)$Symbol <- rownames(sce)
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

    if (file.exists(probabilities)) {
        res <- readRDS(probabilities)
    } else {
    if (length(unique(sce_markers$sample)) > 1) {
        design <- data.frame(id=colnames(sce_markers),sample=sce_markers$sample)
        design <- as.data.table(design)
        design <- model.matrix(~sample, data = design)
        res <- cellassign(exprs_obj = sce_markers, X = design, s = s, marker_gene_info = rho, B = 20, shrinkage = TRUE, verbose = TRUE, rel_tol_em = 1e-5, num_runs=1)
    } else {
        res <- cellassign(exprs_obj = sce_markers, s = s, marker_gene_info = rho, B = 20, shrinkage = TRUE, verbose = TRUE, rel_tol_em = 1e-5, num_runs=1)
    }
    saveRDS(res, file=probabilities)
    }
    sce_markers$cell_type <- res$cell_type
    sce$cell_type <- "Other"
    sce[,colnames(sce_markers)]$cell_type <- sce_markers$cell_type
    merged$cell_type <- sce$cell_type
}

celltypes <- data.frame(barcode=colnames(sce), celltype=sce$cell_type)
write.csv(celltypes, file=celltype_csv)
