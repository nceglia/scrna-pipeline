library(Seurat)
library(SingleCellExperiment)
library(scran)
library(reticulate)
library(cellassign)
library(tensorflow)
library(mltools)
library(data.table)
library(cowplot)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
merged_file      <- args[1]
marker_csv       <- args[2]
probabilities    <- args[3]
annotated_seurat <- args[4]
batch_figure     <- args[5]

merged <- readRDS(merged_file)

DefaultAssay(object = merged) <- "RNA"
sce <- as.SingleCellExperiment(merged)
counts(sce) <- data.matrix(counts(sce))
sce <- computeSumFactors(sce)
sce$size_factor <- sizeFactors(sce)

rho <- read.csv(marker_csv)
rownames(rho) <- rho$X
rho <- rho[,-1]

sce_markers <- sce[rowSums(counts(sce)) > 1000,]
genes <- intersect(rownames(sce_markers),rownames(rho))
sce_markers <- sce_markers[genes,]
sce_markers <- sce_markers[,colSums(counts(sce_markers))>1000]
genes <- intersect(rownames(sce_markers),rownames(rho))
rho <- rho[genes,]
s <- sizeFactors(sce_markers)
rho <- data.matrix(rho)
design <- data.frame(id=colnames(sce_markers),sample=sce_markers$sample)
design <- as.data.table(design)
design <- model.matrix(~sample, data = design)

res <- cellassign(exprs_obj = sce_markers, X = design, s = s, marker_gene_info = rho, B = 20, shrinkage = TRUE, verbose = TRUE, rel_tol_em = 1e-1, num_runs=1)
saveRDS(res, file=probabilities)

sce$cell_type <- "Other"
sce[,colnames(sce_markers)]$cell_type <- res$cell_type
merged$cell_type <- sce$cell_type

merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged, vars.to.regress = c("CC.Diff", "percent.mt", "nCount_RNA"))
merged <- RunPCA(merged, verbose = TRUE)
merged <- RunUMAP(merged, dims = 1:50)

saveRDS(merged, file=annotated_seurat)

Idents(merged) <- "cell_type"
ctumap <- DimPlot(merged,reduction="umap")
Idents(merged) <- "sample"
sumap <- DimPlot(merged,reduction="umap")

figure <- plot_grid(ctumap, sumap, nrow = 1, align = "h")
ggsave(batch_figure,figure,width=8,height=3)
