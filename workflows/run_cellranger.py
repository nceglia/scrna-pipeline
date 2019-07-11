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

def RunCellranger(sampleid, finished, reference):
    CellRanger.count([sampleid],reference_override=reference)
    open(finished,"w").write("Completed")

def RunUpload(sampleid,  finished):
    tenx_output = os.path.join(config.jobpath,"{}/outs/".format(sampleid))
    tenx = TenxAnalysis(tenx_output)
    tenx.finalize()
    tenxds = TenxDataStorage(sampleid)
    tenxds.upload_cellranger(tenx)
    open(finished,"w").write("Completed")


def RunCellranger(sampleid, workflow, full=False):
    workflow.transform (
        name = "cellranger_counts_human",
        func = RunCellranger,
        args = (
            sampleid,
            pypeliner.managed.OutputFile("cellranger_human.complete"),
            config.reference
        )
    )
    workflow.transform (
        name = "cellranger_counts_mouse",
        func = RunCellranger,
        args = (
            sampleid,
            pypeliner.managed.OutputFile("cellranger_mouse.complete"),
            config.mouse_reference
        )
    )
    workflow.transform (
        name = "cellranger_upload_human",
        func = RunUpload,
        args = (
            sampleid,
            pypeliner.managed.InputFile("cellranger_human.complete"),
            pypeliner.managed.OutputFile("human_upload.complete"),
        )
    )
    workflow.transform (
        name = "cellranger_upload_mouse",
        func = RunUpload,
        args = (
            sampleid,
            pypeliner.managed.InputFile("cellranger_mouse.complete"),
            pypeliner.managed.OutputFile("mouse_upload.complete"),
        )
    )
    return workflow
