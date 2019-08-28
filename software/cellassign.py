import subprocess
import numpy
import os
import sys
import pickle
import shutil

"""SINGULARITYENV_NPY_MKL_FORCE_INTEL=GNU SINGULARITYENV_PATH=$PATH:/common/juno/OS7/10.1/linux3.10-glibc2.17-x86_64/bin/ SINGULARITYENV_R_HOME=/usr/local/lib/R/ SINGULARITYENV_LD_LIBRARY_PATH=/usr/local/lib/R/lib"""

from interface.genemarkermatrix import GeneMarkerMatrix

class CellAssign(object):

    @staticmethod
    def cmd(rdata, rho_csv, results, lsf=True):
        CellAssign.script(rdata, rho_csv, results)
        env = os.environ.copy()
        submit = ["Rscript","{}/run_cellassign.R".format(os.path.split(rdata)[0])]
        if lsf:
            submit = """bsub -K -J cellassign -R "rusage[mem=4]" -R "select[type==CentOS7]" -n 16 -We 40 -o out -e err Rscript {}""".format("{}/run_cellassign.R".format(os.path.split(rdata)[0])).split()
        else:
            submit = ["Rscript","{}/run_cellassign.R".format(os.path.split(rdata)[0])]
        subprocess.call(submit, env=env)
        matched_results = os.path.join(os.path.split(results)[0],"cell_types.tsv")
        if lsf:
            submit = """bsub -K -J mapcells -R "rusage[mem=4]" -R "select[type==CentOS7]" -n 16 -We 40 -o out -e err Rscript {}""".format("{}/match.R".format(os.path.split(rdata)[0])).split()
        else:
            submit = ["Rscript","{}/match.R".format(os.path.split(rdata)[0])]
        subprocess.call(submit, env=env)

    @staticmethod
    def run(rdata, rho_yaml, results, rho_csv=".cache/rho.csv"):
        if not os.path.exists(".cache"):
            os.makedirs(".cache")
        marker_list = GeneMarkerMatrix.read_yaml(rho_yaml)
        marker_list.write_matrix(rho_csv)
        assert os.path.exists(rho_csv)
        CellAssign.cmd(rdata, rho_csv, results)
        print ("CellAssign finished.")
        matched_results = os.path.join(os.path.split(rdata)[0],"cell_types.tsv")
        pkl_fit = os.path.join(os.path.split(rdata)[0],"cell_types.pkl")
        lines = open(matched_results,"r").read().splitlines()
        header = lines.pop(0)
        barcodes = []
        celltypes = []
        pyfit = dict()
        for line in lines:
            line = [x.replace('"','') for x in line.split(",")]
            barcodes.append(line[1])
            celltypes.append(line[2])
        pyfit["Barcode"] = barcodes
        pyfit["cell_type"] = celltypes
        pickle.dump(pyfit, open(pkl_fit,"wb"))
        print ("Results written.")

    @staticmethod
    def script(rdata, rho_csv, results):
        filtered_sce = os.path.join(os.path.split(rdata)[0],"sce_cas.rdata")
        filtered_rho = os.path.join(os.path.split(rdata)[0],"rho_cas.rdata")
        matched_results = os.path.join(os.path.split(results)[0],"cell_types.tsv")
        configured = open("{}/run_cellassign.R".format(os.path.split(rdata)[0]),"w")
        configured.write(script.format(sce=rdata,rho=rho_csv,fname=results,fsce=filtered_sce,frho=filtered_rho))
        configured.close()
        match = open("{}/match.R".format(os.path.split(rdata)[0]),"w")
        match.write(match_barcodes.format(sce=filtered_sce,fit=results,fname=matched_results))
        match.close()


script = """
library(reticulate)
use_python("/home/ceglian/anaconda/bin/python3")
library(cellassign)
library(tensorflow)
library(scran)

rho <- read.csv("{rho}")
rownames(rho) <- rho$X
rho <- rho[,-1]

sce <- readRDS("{sce}")
colnames(sce) <- sce$Barcode
sce_result <- sce

print('filter')
rownames(sce) <- rowData(sce)$Symbol
rho <- as.matrix(rho)
counts(sce) <- data.matrix(counts(sce))
sce <- sce[rowSums(counts(sce)) > 0,]
common_genes <- intersect(rowData(sce)$Symbol,rownames(rho))
sce <- sce[common_genes,]
sce <- sce[,colSums(counts(sce))>0]
common_genes <- intersect(rowData(sce)$Symbol,rownames(rho))
rho <- rho[common_genes,]


rho <- data.matrix(rho)
s <- sizeFactors(sce)
print('call')
fit_cellassign <- cellassign(exprs_obj = sce, marker_gene_info = rho, s = s, B=20, shrinkage=TRUE)
sce_result <- sce_result[,colnames(sce)]
colData(sce_result)$cell_type <- fit_cellassign$cell_type
print('save')
saveRDS(fit_cellassign, file = '{fname}')
saveRDS(sce_result, file="{fsce}")
saveRDS(rho, file="{frho}")
"""


match_barcodes = """
library(SingleCellExperiment)
library(cellassign)
sce = readRDS("{sce}")
fit = readRDS("{fit}")

barcodes <- data.frame(barcode=colData(sce)$Barcode,celltype=fit$cell_type)
write.csv(barcodes, file="{fname}")

"""
