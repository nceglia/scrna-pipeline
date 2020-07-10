library(Seurat)
library(SingleCellExperiment)
library(DropletUtils)

args = commandArgs(trailingOnly=TRUE)
matrix <- args[1]
mito   <- args[2]
object_file <- args[3]
sce_file <- args[4]
sce_raw  <- args[5]
out_matrix <- args[6]

if (grepl(".rdata", matrix, fixed = TRUE)) {
    sce <- readRDS(matrix)
    matrix <- counts(sce)
    colnames(matrix) <- colData(sce)$Barcode
    rownames(matrix) <- rowData(sce)$Symbol
    write10xCounts(out_matrix, matrix)
    data <- Read10X(data.dir = out_matrix)
} else {
    data <- Read10X(data.dir = matrix)
}

seurat <- CreateSeuratObject(counts = data, project = "RNASCP", min.cells = 3, min.features = 500)
sce <- as.SingleCellExperiment(seurat)
saveRDS(sce, file=sce_raw)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, pattern = "^RPS|^RPL")
seurat <- subset(seurat, subset = nFeature_RNA > 500 & percent.mt < mito & nCount_RNA > 1000)
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
saveRDS(seurat, file=object_file)
sce <- as.SingleCellExperiment(seurat)
saveRDS(sce, file=sce_file)
