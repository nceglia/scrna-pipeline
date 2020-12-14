library(Seurat)
library(SingleCellExperiment)
library(DropletUtils)
library(scater)

args = commandArgs(trailingOnly=TRUE)
matrix <- args[1]
mito   <- args[2]
ncounts <- args[3]
object_file <- args[4]
sce_file <- args[5]
sce_raw  <- args[6]
out_matrix <- args[7]
nfeatures <- args[8]

if (grepl(".rdata", matrix, fixed = TRUE)) {
    sce <- readRDS(matrix)
    rownames(sce) <- rowData(sce)$Symbol
    symbols <- lapply(rownames(sce), function(x) { toupper(x) })
    rownames(sce) <- symbols
    rowData(sce)$Symbol <- rownames(sce)
    # seurat <- CreateSeuratObject(counts = data, project = "RNASCP", min.cells = 3)
    seurat <- CreateSeuratObject(counts(sce), project = "RNASCP", min.cells=3)
    # matrix <- counts(sce)
    # colnames(matrix) <- colData(sce)$Barcode
    # rownames(matrix) <- rowData(sce)$Symbol
    # write10xCounts(out_matrix, matrix)
    # data <- Read10X(data.dir = out_matrix)
} else if (grepl(".txt", matrix, fixed = TRUE)) {
    x <- read.table(args[1])
    seurat <- CreateSeuratObject(x, project="RNASCP", min.cells=3)
    sce <- as.SingleCellExperiment(seurat)
} else {
    print(matrix)
    data <- Read10X(data.dir = matrix)
    seurat <- CreateSeuratObject(counts = data, project = "RNASCP", min.cells = 3)
    sce <- as.SingleCellExperiment(seurat)
}
saveRDS(sce, file=sce_raw)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, pattern = "^RPS|^RPL")
seurat <- subset(seurat, subset = nFeature_RNA > 100 & percent.mt <= 30 & nCount_RNA > 100)
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
saveRDS(seurat, file=object_file)
sce <- as.SingleCellExperiment(seurat)
saveRDS(sce, file=sce_file)
