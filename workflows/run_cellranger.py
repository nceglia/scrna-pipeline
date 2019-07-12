import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os
import glob

from software.cellranger import CellRanger
from interface.fastqdirectory import FastQDirectory
from interface.tenxanalysis import TenxAnalysis
from utils.cloud import TenxDataStorage
from utils.config import Configuration, write_config

config = Configuration()

def DownloadFastqs(sampleid, finished):
    fastqs = [FastQDirectory(fastq, config.prefix, config.jobpath, config.datapath) for fastq in [sampleid]]
    fastqs = glob.glob(os.path.join(config.datapath,"*.fastq.gz"))
    assert len(fastqs) > 0, "No Fastqs Download or Found."
    open(finished,"w").write("Completed")


def Counts(sampleid, fastqs_downloaded, finished, reference):
    CellRanger.count([sampleid],reference_override=reference)
    open(finished,"w").write("Completed")

def RunUpload(sampleid, before, finished, species):
    if not os.path.exists(finished):
        tenx_output = os.path.join(config.jobpath,"{}/outs/".format(sampleid))
        tenx = TenxAnalysis(tenx_output)
        tenx.finalize()
        tenxds = TenxDataStorage(sampleid, species=species)
        tenxds.upload_cellranger(tenx)
    open(finished,"w").write("Completed")

def RunCellranger(sampleid, workflow):
    workflow.transform (
        name = "download_fastqs",
        func = DownloadFastqs,
        args = (
            sampleid,
            pypeliner.managed.TempOutputFile("download_fastqs.complete"),
        )
    )
    workflow.transform (
        name = "cellranger_counts_human",
        func = Counts,
        args = (
            sampleid,
            pypeliner.managed.TempInputFile("download_fastqs.complete"),
            pypeliner.managed.TempOutputFile("cellranger_human.complete"),
            config.reference
        )
    )
    workflow.transform (
        name = "cellranger_counts_mouse",
        func = Counts,
        args = (
            sampleid,
            pypeliner.managed.TempInputFile("download_fastqs.complete"),
            pypeliner.managed.TempOutputFile("cellranger_mouse.complete"),
            config.mouse_reference
        )
    )
    workflow.transform (
        name = "cellranger_upload_human",
        func = RunUpload,
        args = (
            sampleid,
            pypeliner.managed.TempInputFile("cellranger_human.complete"),
            pypeliner.managed.TempOutputFile("human_upload.complete"),
            "human"
        )
    )
    workflow.transform (
        name = "cellranger_upload_mouse",
        func = RunUpload,
        args = (
            sampleid + "_mouse",
            pypeliner.managed.TempInputFile("cellranger_mouse.complete"),
            pypeliner.managed.TempOutputFile("mouse_upload.complete"),
            "mouse"
        )
    )
    return workflow
