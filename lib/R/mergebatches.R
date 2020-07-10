library(SingleCellExperiment)
library(Seurat)

args = commandArgs(trailingOnly=TRUE)
paths <- read.csv(args[1], sep="\t", header=TRUE)
object_file <- args[2]
celltype_csv <- args[3]

i <- 0
rdata <- as.character(paths$path)
batches <- as.character(paths$batch_id)
seurats <- vector("list", length(rdata)-1)

init_seurat <- readRDS(rdata[[1]])
init_seurat$batch <-batches[1]
for (i in 2:length(rdata)) {
  seurat <- readRDS(rdata[i])
  seurat$batch <- batches[i]
  seurats[[i-1]] <- seurat
}
merged <- merge(init_seurat, y = seurats, project = "RNASCP")

merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged, vars.to.regress = c("CC.Diff", "percent.mt", "nCount_RNA"))
merged <- RunPCA(merged, verbose = TRUE)
merged <- RunUMAP(merged, dims = 1:50)

saveRDS(merged, file=object_file)
# merged <- readRDS(object_file)

write.csv(data.frame(table(merged$cell_type)), file=celltype_csv)
