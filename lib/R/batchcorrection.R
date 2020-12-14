library(Seurat)
library(harmony)
library(scater)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(schex)

options(future.globals.maxSize= 2097152000)

args = commandArgs(trailingOnly=TRUE)

seurat <- ReadH5AD(args[1])
integrated_object_file <- args[2]
project_figure <- args[3]
celltype_csv <- args[4]
cell_annotations <- args[5]
metadata <- args[6]

df <- read.csv(metadata)

rownames(df) <- df$cell_id

seurat <- AddMetaData(seurat, df$cell_type, col.name = "cell_type")
seurat <- AddMetaData(seurat, df$super_cell_type, col.name = "super_cell_type")
seurat <- AddMetaData(seurat, df$nCount_RNA, col.name = "nCount_RNA")
seurat <- AddMetaData(seurat, df$nFeature_RNA, col.name = "nFeature_RNA")
seurat <- AddMetaData(seurat, df$percent.mt, col.name = "percent.mt")
seurat <- AddMetaData(seurat, df$percent.rb, col.name = "percent.rb")
seurat <- AddMetaData(seurat, df$doublet, col.name = "doublet")
seurat <- AddMetaData(seurat, df$doublet_score, col.name = "doublet_score")
seurat <- AddMetaData(seurat, df$sample, col.name = "sample")
seurat <- AddMetaData(seurat, df$Phase, col.name = "Phase")
seurat <- AddMetaData(seurat, df$batch, col.name = "batch")


merged <- NormalizeData(seurat)
merged <- FindVariableFeatures(merged, nfeatures = 3000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, verbose = F)
RunUMAP_wrapper <- function(seurat_object, ndim) RunUMAP(seurat_object, dims = 1:ndim, reduction.name = paste0("umap", ndim), reduction.key = paste0("umap", ndim, "_"))
merged <- RunUMAP_wrapper(merged, 50)
saveRDS(merged, file="/work/shah/ceglian/cohort_all.rds")
merged <- FindNeighbors(merged, dims = 1:50)

cell_types_super <- c(T.cell = "T.super", B.cell = "B.super", Plasma.cell = "B.super", Mast.cell = "Myeloid.super", Dendritic.cell = "Myeloid.super", Monocyte = "Myeloid.super", Fibroblast = "Stromal.super", Endothelial.cell = "Stromal.super", Ovarian.cancer.cell = "Ovarian.cancer.super", Other = "Other.super")
merged$cell_type_super <- cell_types_super[merged$cell_type]

# ### cohort outputs
write_rds(merged, integrated_object_file)
seu_slim <- as_tibble(cbind(cell_id = colnames(merged), FetchData(merged, c("cell_type", "cell_type_super", "batch", "percent.mt", "nCount_RNA", "nFeature_RNA", "doublet_score", "doublet", "sample","umap50_1", "umap50_2"))))
write_tsv(seu_slim, cell_annotations)