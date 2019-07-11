import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os
import shutil

from interface.tenxanalysis import TenxAnalysis
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl

from utils.config import Configuration, write_config

config = Configuration()

def Run(sampleid, species, umi_plot, mito_plot, ribo_plot, counts_plot, raw_sce):
    print("Running QC.")
    tenx = TenxDataStorage(sampleid, version="v3", species=species)
    tenx.download()
    tenx_analysis = TenxAnalysis(tenx.tenx_path)
    tenx_analysis.load()
    tenx_analysis.extract()
    print("Extracted.")
    qc = QualityControl(tenx_analysis,sampleid)
    if not os.path.exists(qc.sce):
        qc.run(mito=config.mito)
        print ("Uploading")
        qc.upload_raw()
        qc.upload()
    plots = qc.plots
    umi = os.path.join(plots,"umi.png")
    mito = os.path.join(plots,"mito.png")
    ribo = os.path.join(plots, "ribo.png")
    counts = os.path.join(plots, "counts.png")
    cvf = os.path.join(plots, "total_counts_v_features.png")

    results = os.path.join(config.jobpath, "results")
    if not os.path.exists(results):
        os.makedirs(results)

    shutil.copyfile(umi, umi_plot)
    shutil.copyfile(mito, mito_plot)
    shutil.copyfile(ribo, ribo_plot)
    shutil.copyfile(counts, counts_plot)
    shutil.copyfile(qc.sce, raw_sce)

def RunQC(sampleid, workflow, species=None):
    if species == "mouse":
        umi = "mouse_umi.png"
        mito = "mouse_mito.png"
        ribo = "mouse_ribo.png"
        counts = "mouse_counts.png"
        sce = "mouse_raw_sce.rdata"
    else:
        umi = "umi.png"
        mito = "mito.png"
        ribo = "ribo.png"
        counts = "counts.png"
        sce = "raw_sce.rdata"
    workflow.transform (
        name = "quality_control_{}".format(species),
        func = Run,
        args = (
            sampleid,
            species,
            # pypeliner.managed.TempInputFile("cellranger_{}.complete".format(species)),
            pypeliner.managed.TempOutputFile(umi),
            pypeliner.managed.TempOutputFile(mito),
            pypeliner.managed.TempOutputFile(ribo),
            pypeliner.managed.TempOutputFile(counts),
            pypeliner.managed.TempOutputFile(sce),
        )
    )
    return workflow
