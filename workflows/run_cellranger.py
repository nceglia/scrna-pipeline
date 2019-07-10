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

def RunCellranger(sampleid, finished, ref_override, full):
    if full:
        CellRanger.count([sampleid],ref_override=ref_override)
    open(finished,"w").write("Completed")

def RunUpload(sampleid, finished):
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
            False,
            full
        )
    )
    workflow.transform (
        name = "cellranger_counts_mixed_ref",
        func = RunCellranger,
        args = (
            sampleid,
            pypeliner.managed.OutputFile("cellranger_mixed_ref.complete"),
            True,
            full
        )
    )
    workflow.transform (
        name = "cellranger_upload",
        func = RunUpload,
        args = (
            sampleid,
            pypeliner.managed.InputFile("cellranger.complete"),
            pypeliner.managed.InputFile("cellranger_counts_mixed_ref.complete"),
            pypeliner.managed.OutputFile("upload.complete"),
        )
    )
    return workflow
