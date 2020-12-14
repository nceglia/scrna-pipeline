from singlecellexperiment import SingleCellExperiment
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

sce = SingleCellExperiment.fromRData("pdx_filtered_sce.rdata")
genes = sce.rowData["Symbol"]
print(sce.assays["counts"][0])

# gene_index = dict(zip(genes,list(range(len(genes))))) #{"VIM": 0, "PTPRC": 1}

# embeds = nn.Embedding(len(genes), 30)
