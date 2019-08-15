import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os
import json
import shutil
import subprocess
from collections import defaultdict
from scipy import stats
import math
import numpy

from interface.tenxanalysis import TenxAnalysis
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl
from utils.cloud import SampleCollection
from interface.qualitycontrol import QualityControl
from interface.genemarkermatrix import GeneMarkerMatrix
from utils.plotting import celltypes, tsne_by_cell_type, umap_by_cell_type
from software.cellassign import CellAssign

from utils.config import Configuration, write_config

config = Configuration()

def RunDownload(sampleids, finished):
    for i, sample in enumerate(sampleids):
        tenx = TenxDataStorage(sample)
        path = tenx.download()
        path_json = {sample: path}
        open(finished(i),"w").write(json.dumps(path_json))

def RunExtract(sample_to_path, rdata_path, summary_path):
    sample = json.loads(open(sample_to_path,"r").read())
    sampleid, path = list(sample.items()).pop()
    tenx_analysis = TenxAnalysis(path)
    tenx_analysis.load()
    tenx_analysis.extract()
    qc = QualityControl(tenx_analysis, sampleid)
    if not os.path.exists(qc.sce):
        qc.run(mito=config.mito)
    shutil.copyfile(tenx_analysis.summary, summary_path)
    shutil.copyfile(qc.sce, rdata_path)

def RunCellAssign(sce, annot_sce):
    _rho_csv = os.path.join(os.path.split(sce)[0],"rho_csv_sub.csv")
    _fit = os.path.join(os.path.split(sce)[0],"fit_sub.pkl")
    sampleid = sce.split("/")[-2]
    #CellAssign.run(sce, config.rho_matrix, _fit, rho_csv=_rho_csv)
    filtered_sce = os.path.join(os.path.split(sce)[0],"sce_cas.rdata")
    shutil.copyfile(filtered_sce, annot_sce)

def RunModeCopyNumber(copy_number_data):
    output = open(copy_number_data,"w")
    cell_to_clone = open(config.copy_cell_clones,"r").read().splitlines()
    cell_to_clone.pop(0)
    cell_to_clone = dict([x.split("\t") for x in cell_to_clone])
    copy_number_data = open(config.filtered_cell_cn,"r").read().splitlines()
    header = copy_number_data.pop(0).split("\t")
    print("chr,start,end,copy_number,clone")
    copy_number_mapping = defaultdict(lambda : defaultdict(list))
    for row in copy_number_data:
            row = row.split("\t")
            chrm, start, stop, width = row[:4]
            stop = str(int(start) + int(width))
            cells = dict(zip(header[4:],row[4:]))
            #ploidy = stats.mode([math.ceil(int(cn)) for cn in cells.values()])
            #ploidy = ploidy.mode[0]
            #if ploidy == 0:
            #       ploidy = 2
            for cell, copy_number in cells.items():
                #normalized_cn = math.ceil(int(copy_number) / (int(ploidy) / 2))
                copy_number_mapping[chrm+":"+start+":"+stop][cell_to_clone[cell]].append(int(copy_number))

    for loc in copy_number_mapping:
            chrm, start, end = loc.split(":")
            for clone, copy_numbers in copy_number_mapping[loc].items():
                    if clone == "None": continue
                    mode_copy_number = stats.mode(copy_numbers).mode[0]
                    line = ["chr"+chrm, start, end, str(int(mode_copy_number)), clone]
                    print(",".join(line))
                    output.write(",".join(line)+"\n")
    output.close()

def RunCloneAlignInput(sce, copy_number_data, clone_sce, cnv_mat):
    seurat_cached = os.path.join("seurat_raw.rdata")
    sce_cached = os.path.join(os.path.split(sce)[0],"sce_cas.rdata")
    rcode = """
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(dplyr)
    library(tidyverse)
    library(org.Hs.eg.db)
    library(SingleCellExperiment)

    sce <- readRDS('{sce}')
    sce <- sce[,sce$cell_type=="Ovarian.cancer.cell"]
    rownames(sce) <- rowData(sce)$ensembl_gene_id
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    g <- genes(txdb, single.strand.genes.only=FALSE)
    entrezgene_ensembl_map <- as.list(org.Hs.egENSEMBL)
    entrezgene_ensembl_map <- lapply(entrezgene_ensembl_map, `[`, 1)
    cnv_data <- read_csv(file='{cnv}')
    cnv_gr <- makeGRangesFromDataFrame(cnv_data, keep.extra.columns = TRUE)
    olaps <- findOverlaps(g, cnv_gr)
    df_gene <- data_frame(entrezgene = names(g)[queryHits(olaps)],
               copy_number = mcols(cnv_gr)$copy_number[subjectHits(olaps)],
               clone = mcols(cnv_gr)$clone[subjectHits(olaps)])
    df_gene <- dplyr::filter(df_gene, entrezgene %in% names(entrezgene_ensembl_map)) %>%
      dplyr::mutate(ensembl_gene_id = unlist(entrezgene_ensembl_map[entrezgene])) %>%
      dplyr::select(ensembl_gene_id, entrezgene, copy_number, clone) %>%
      drop_na()
    df_gene <- df_gene %>% group_by(ensembl_gene_id) %>% tally() %>%
      filter(n == length(unique(df_gene$clone))) %>%
      inner_join(df_gene) %>%
      dplyr::select(-n)
    df_gene_expanded <- spread(df_gene, clone, copy_number)
    cnv_mat <- dplyr::select(df_gene_expanded, -ensembl_gene_id, -entrezgene) %>% as.matrix()
    rownames(cnv_mat) <- df_gene_expanded$ensembl_gene_id
    keep_gene <- rowMins(cnv_mat) <= 6 & rowVars(cnv_mat) > 0
    cnv_mat <- cnv_mat[keep_gene,]
    common_genes <- intersect(rownames(cnv_mat), rownames(sce))
    sce <- sce[common_genes,]
    cnv_mat <- cnv_mat[common_genes,]
    saveRDS(file='{clone_sce}', sce)
    saveRDS(file='{cnv_mat}', cnv_mat)
    """
    path = os.path.split(sce)[0]
    convert_script = os.path.join(path,"build_input.R")
    output = open(convert_script,"w")
    output.write(rcode.format(sce=sce,cnv=copy_number_data,clone_sce=clone_sce,cnv_mat=cnv_mat))
    output.close()
    subprocess.call(["Rscript","{}".format(convert_script)])

def RunCloneAlign(clone_sce, cnv_mat, annotated_sce, cal_fit):
    rdata = """
    library(clonealign)
    sce <- readRDS('{clone_sce}')
    cnv <- readRDS('{cnv_mat}')
    cal <- clonealign(sce, cnv)
    sce$clone <- cal$clone
    saveRDS(sce,file='{annotated_sce}')
    saveRDS(cal, file='{cal_fit}')
    """
    path = os.path.split(clone_sce)[0]
    run_script = os.path.join(path,"run_clonealign.R")
    output.write(rcode.format(clone_sce=clone_sce,cnv_mat=cnv_mat,annotated_sce=annotated_sce,cal_fit=cal_fit))
    output.close()
    subprocess.call(["Rscript","{}".format(run_script)])


def RunCloneAlignWorkflow(workflow):
    print("Creating workflow.")
    all_samples = [config.prefix]
    workflow.transform (
        name = "download_collection",
        func = RunDownload,
        args = (
            all_samples,
            pypeliner.managed.TempOutputFile("sample_path.json","sample")
        )
    )
    workflow.transform (
        name = "extract_rdata",
        func = RunExtract,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("sample_path.json","sample"),
            pypeliner.managed.TempOutputFile("sample.rdata","sample"),
            pypeliner.managed.TempOutputFile("summary_path.html","sample")
        )
    )
    workflow.transform (
        name = "run_cellassign",
        func = RunCellAssign,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("sample.rdata","sample"),
            pypeliner.managed.TempOutputFile("sce.rdata","sample"),
        )
    )
    workflow.transform (
        name = "run_modecn",
        func = RunModeCopyNumber,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempOutputFile("copy_number_data.csv","sample"),
        )
    )
    workflow.transform (
        name = "run_clonealigninputs",
        func = RunCloneAlignInput,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("sce.rdata","sample"),
            pypeliner.managed.TempInputFile("copy_number_data.csv","sample"),
            pypeliner.managed.TempOutputFile("clone_sce.rdata","sample"),
            pypeliner.managed.TempOutputFile("cnv.rdata","sample"),
        )
    )
    workflow.transform (
        name = "run_clonealign",
        func = RunCloneAlign,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("clone_sce.rdata","sample"),
            pypeliner.managed.TempInputFile("cnv.rdata","sample"),
            pypeliner.managed.TempOutputFile("annotated.rdata","sample"),
            pypeliner.managed.TempOutputFile("cal.rdata","sample"),
        )
    )
    # workflow.transform (
    #     name = "run_figures",
    #     func = RunCloneFigures,
    #     axes = ('sample',),
    #     args = (
    #         pypeliner.managed.TempInputFile("sample.rdata","sample"),
    #         pypeliner.managed.TempOutputFile("sce.rdata","sample"),
    #     )
    # )
    return workflow
