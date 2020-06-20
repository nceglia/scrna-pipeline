library(Seurat)
library(SingleCellExperiment)

args = commandArgs(trailingOnly=TRUE)
filtered_matrices <- args[1]
object_file <- args[2]

data <- Read10X(data.dir = filtered_matrices)
seurat <- CreateSeuratObject(counts = data, project = "SPECTRUM", min.cells = 3, min.features = 500)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, pattern = "^RPS|^RPL")
seurat <- subset(seurat, subset = nFeature_RNA > 500 & percent.mt < 25 & nCount_RNA > 1000)
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
saveRDS(seurat, file=object_file)

