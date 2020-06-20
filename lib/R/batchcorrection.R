library(Seurat)
library(harmony)
library(scater)
library(tidyverse)

options(future.globals.maxSize= 2097152000)

args = commandArgs(trailingOnly=TRUE)

merged <- readRDS(args[1])
integrated_object_file <- args[2]

integrated <- RunHarmony(merged, group.by.vars="batch", dims.use = 1:50)
integrated <- RunPCA(integrated, verbose = TRUE)
integrated <- RunUMAP(integrated, dims = 1:50, reduction = "harmony", reduction.name = "umapharmony", reduction.key = "umapharmony_")
saveRDS(integrated, file=integrated_object_file)
