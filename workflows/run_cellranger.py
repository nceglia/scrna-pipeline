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

def RunCellranger(sampleid, finished, full):
    if full:
        CellRanger.count([sampleid])
    open(finished,"w").write("Completed")

def RunUpload(sampleid, run, finished):
    tenx_output = os.path.join(config.jobpath,"{}/outs/".format(sampleid))
    tenx = TenxAnalysis(tenx_output)
    tenx.finalize()
    tenxds = TenxDataStorage(sampleid)
    tenxds.upload_cellranger(tenx)
    open(finished,"w").write("Completed")


def RunCellranger(sampleid, workflow, full=False):
    workflow.transform (
        name = "cellranger_counts",
        func = RunCellranger,
        args = (
            sampleid,
            pypeliner.managed.OutputFile("cellranger.complete"),
            full
        )
    )
    workflow.transform (
        name = "cellranger_upload",
        func = RunUpload,
        args = (
            sampleid,
            pypeliner.managed.InputFile("cellranger.complete"),
            pypeliner.managed.OutputFile("upload.complete"),
        )
    )
    return workflow
