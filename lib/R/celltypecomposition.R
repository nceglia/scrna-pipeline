library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
seurat <- readRDS(args[1])
batch_composition <- args[2]
sample_composition <- args[3]

df <- data.frame(cell_type=seurat$cell_type, sample=seurat$sample)
sample_summary <- summarise(group_by(df,cell_type,sample),count =n())
sample_summary <- data.frame(sample_summary)

df <- data.frame(cell_type=seurat$cell_type, batch=seurat$batch)
batch_summary <- summarise(group_by(df,cell_type,batch),count =n())
batch_summary <- data.frame(batch_summary)


fig <- ggplot(sample_summary, aes(factor(sample), count, fill = cell_type)) + geom_bar(stat="identity",position="fill") + xlab("Sample") + ylab("Cell Type") + theme(legend.title = element_blank()) + ggtitle("")
ggsave(sample_composition, fig,width=4,height=3)


fig <- ggplot(batch_summary, aes(factor(batch), count, fill = cell_type)) + geom_bar(stat="identity",position="fill") + xlab("Batch") + ylab("Cell Type") + theme(legend.title = element_blank()) + ggtitle("")
ggsave(batch_composition, fig,width=4,height=3)
