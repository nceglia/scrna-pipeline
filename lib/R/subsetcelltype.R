library(Seurat)
library(scater)
library(tidyverse)
library(harmony)
library(cowplot)
library(ggplot2)
library(schex)
library(SingleCellExperiment)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)    
merged <- readRDS(args[1])
celltype <- args[2]
celltype_obj <- args[3]
marker_csv <- args[4]
markers_tsv <- args[5]
cells_tsv <- args[6]
resolution <- as.double(args[7])

message("BLAHJaasdfadsgadg")
message(resolution)
message("BLAFADSFSDS")
Idents(merged) <- "super_cell_type"
seu_obj <- merged[,merged$super_cell_type==celltype]
if (length(colnames(seu_obj)) < 100) {
  quit()
}

cc.genes <- read_lines("/work/shah/ceglian/isabl2/shahlab_apps/shahlab_apps/apps/cellassign/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43] ## s-phase genes
g2m.genes <- cc.genes[44:97] ## g2m transition genes
score_cc <- function(sobj) {
sobj <- CellCycleScoring(sobj, s.genes, g2m.genes)
sobj$CC.Diff <- sobj$S.Score - sobj$G2M.Score
return(sobj)
}

seu_obj$patient_id <- str_sub(seu_obj$sample, 0, 15)
seu_obj <- NormalizeData(seu_obj)
seu_obj <- FindVariableFeatures(seu_obj)
seu_obj <- ScaleData(seu_obj)
seu_obj <- score_cc(seu_obj)
seu_obj <- ScaleData(seu_obj, vars.to.regress = c("CC.Diff", "percent.mt", "nCount_RNA"))
seu_obj <- RunPCA(seu_obj, verbose = F)
seu_obj <- RunUMAP(seu_obj, dims = 1:50)
seu_obj <- RunHarmony(seu_obj, group.by.vars = "patient_id", dims.use = 1:50, assay.use = "RNA")
seu_obj <- RunUMAP(seu_obj, dims = 1:50, reduction = "harmony", reduction.name = "umapharmony", reduction.key = "umapharmony_")
seu_obj <- FindNeighbors(seu_obj, dims = 1:50, reduction = "harmony")
seu_obj <- FindClusters(seu_obj, resolution = resolution)

saveRDS(seu_obj, file=celltype_obj)

sce <- as.SingleCellExperiment(seu_obj)
if ("RNA_snn_res.0.1" %in% colnames(colData(sce))) {
    Idents(seu_obj) <- seu_obj$RNA_snn_res.0.1
    rescol <- "RNA_snn_res.0.1"
} else if ("RNA_snn_res.0.2" %in% colnames(colData(sce))) {
    Idents(seu_obj) <- seu_obj$RNA_snn_res.0.2
    rescol <- "RNA_snn_res.0.2"
} else if ("RNA_snn_res.0.3" %in% colnames(colData(sce))) {
    Idents(seu_obj) <- seu_obj$RNA_snn_res.0.3
    rescol <- "RNA_snn_res.0.3"
}
write_tsv(cbind(tibble(cell_id = colnames(seu_obj)), FetchData(seu_obj, c("sample", "UMAP_1", "UMAP_2", "umapharmony_1", "umapharmony_2", rescol, "cell_type"))), cells_tsv)

seu_markers <- mutate(as_tibble(FindAllMarkers(seu_obj, only.pos = T)), resolution = resolution)
write_tsv(seu_markers, markers_tsv)