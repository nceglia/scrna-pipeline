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

def Run(sampleid, finished):
    tenx_output = os.path.join(config.jobpath,"{}/outs/".format(sampleid))
    tenx = TenxAnalysis(tenx_output)
    tenx.finalize()
    tenxds = TenxDataStorage(sampleid)
    tenxds.upload_cellranger(tenx)
    open(finished,"w").write("Completed")


def RunCellranger(sampleid, workflow):
    workflow.transform (
        name = "cellranger_counts",
        func = Run,
        args = (
            sampleid,
            pypeliner.managed.OutputFile("cellranger.complete"),
        )
    )
    return workflow
