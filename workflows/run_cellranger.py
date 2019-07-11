import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os

from software.cellranger import CellRanger
from interface.fastqdirectory import FastQDirectory
from interface.tenxanalysis import TenxAnalysis
from utils.cloud import TenxDataStorage
from utils.config import Configuration, write_config

config = Configuration()

def Counts(sampleid, finished, reference):
    CellRanger.count([sampleid],reference_override=reference)
    open(finished,"w").write("Completed")

def RunUpload(sampleid, before, finished, species):
    tenx_output = os.path.join(config.jobpath,"{}/outs/".format(sampleid))
    tenx = TenxAnalysis(tenx_output)
    tenx.finalize()
    tenxds = TenxDataStorage(sampleid, species=species)
    tenxds.upload_cellranger(tenx)
    open(finished,"w").write("Completed")

def RunCellranger(sampleid, workflow):
    workflow.transform (
        name = "cellranger_counts_human",
        func = Counts,
        args = (
            sampleid,
            pypeliner.managed.TempOutputFile("cellranger_human.complete"),
            config.reference
        )
    )
    workflow.transform (
        name = "cellranger_counts_mouse",
        func = Counts,
        args = (
            sampleid,
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
