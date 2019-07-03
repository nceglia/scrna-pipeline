import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os
import json
import shutil

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

def RunExtract(sample_to_path, rdata_path):
    sample = json.loads(open(sample_to_path,"r").read())
    sampleid, path = list(sample.items()).pop()
    tenx_analysis = TenxAnalysis(path)
    tenx_analysis.load()
    tenx_analysis.extract()
    qc = QualityControl(tenx_analysis, sampleid)
    if not os.path.exists(qc.sce):
        qc.run(mito=config.mito)
    shutil.copyfile(qc.sce, rdata_path)

def RunCellAssign(sce, annot_sce, rho_csv, fit):
    sampleid = sce.split("/")[-2]
    CellAssign.run(sce, config.rho_matrix, fit, rho_csv=rho_csv)
    filtered_sce = os.path.join(os.path.split(sce)[0],"sce_cas.rdata")
    shutil.copyfile(filtered_sce, annot_sce)

def RunConvert(sce, seurat):
    rcode = """
    library(Seurat)
    library(SingleCellExperiment)
    sce <- readRDS('{sce}')
    seurat <- Convert(from  =sce, to = "seurat")
    saveRDS(seurat,file='{seurat}')
    """
    path = os.path.split(sce)[0]
    convert_script = os.path.join(path,"convert.R")
    output = open(convert_script.format(sce=sce,seurat=seurat),"w")
    output.write(rcode)
    output.close()
    result = subprocess.check_output("R {}".format(convert_script))

def RunSeuratWorkflow(seurat, qcd_seurat):
    rcode = """
    library(Seurat)
    seurat <- NormalizeData(object = {seurat})
    seurat <- FindVariableFeatures(object = seurat)
    seurat <- ScaleData(object = seurat)
    seurat <- RunPCA(object = seurat)
    seurat <- FindNeighbors(object = seurat)
    seurat <- FindClusters(object = seurat)
    seurat <- RunTSNE(object = seurat)
    seurat <- RunUMAP(object = seurat)
    saveRDS(seurat, file = {qcd_seurat})"""
    path = os.path.split(seurat)[0]
    qc_script = os.path.join(path,"qc.R")
    output = open(qc_script.format(seruat=seurat, qcd_seurat=qcd_seurat),"w")
    output.write(rcode)
    output.close()
    result = subprocess.check_output("R {}".format(qc_script))

def RunSeuratViz(seurat, qcd_seurat):
    rcode = """
    library(Seurat)
    DimPlot(object = pbmc, reduction = "tsne")
    DimPlot(object = pbmc, reduction = "pca")
    saveRDS(seurat, file = {qcd_seurat})"""
    path = os.path.split(seurat)[0]
    qc_script = os.path.join(path,"qc.R")
    output = open(qc_script.format(seruat=seurat, qcd_seurat=qcd_seurat),"w")
    output.write(rcode)
    output.close()
    result = subprocess.check_output("R {}".format(qc_script))

def RunCollect(rdata, manifest):
    output = open(manifest,"w")
    for id, rdata in rdata.items():
        output.write(id +"\t" + rdata+"\n")
    output.close()

def RunCollection(workflow):
    all_samples = open(config.samples, "r").read().splitlines()
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
            pypeliner.managed.TempOutputFile("sample.rdata","sample")
        )
    )

    workflow.transform (
        name = "run_cellassign",
        func = RunCellAssign,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("sample.rdata","sample"),
            pypeliner.managed.TempOutputFile("sce.rdata","sample"),
            pypeliner.managed.TempOutputFile("rho.csv","sample"),
            pypeliner.managed.TempOutputFile("cellassignfit.rdata","sample"),
        )
    )

    # workflow.transform (
    #     name = "run_convert",
    #     func = RunConvert,
    #     axes = ('sample',),
    #     args = (
    #         pypeliner.managed.TempInputFile("sce.rdata","sample"),
    #         pypeliner.managed.TempOutputFile("seurat.rdata","sample"),
    #     )
    # )
    #
    # workflow.transform (
    #     name = "run_qc",
    #     func = RunSeuratWorkflow,
    #     axes = ('sample',),
    #     args = (
    #         pypeliner.managed.TempInputFile("seurat.rdata","sample"),
    #         pypeliner.managed.TempOutputFile("seurat_qcd.rdata","sample"),
    #     )
    # )
    #
    # workflow.transform (
    #     name = "visualize_sample",
    #     func = RunSeuratViz,
    #     axes = ('sample',),
    #     args = (
    #         pypeliner.managed.TempInputFile("seurat.rdata","sample"),
    #         pypeliner.managed.TempOutputFile("seurat_qcd.rdata","sample"),
    #     )
    # )

    return workflow
