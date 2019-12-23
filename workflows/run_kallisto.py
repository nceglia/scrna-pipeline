import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os

from software.kallisto import Kallisto
from interface.fastqdirectory import FastQDirectory
from utils.config import Configuration, write_config

config = Configuration()

os.environ["PATH"] = os.environ["PATH"] + ":" + config.kallisto + ":" + config.bustools

def RunKallisto(sampleid, finished):
    fastqs = [FastQDirectory(sampleid, config.prefix, config.jobpath, config.datapath)]
    krunner = Kallisto(fastqs[0], sampleid)
    tenx_path = krunner.count()
    # tenx = TenxAnalysis(tenx_path)
    # tenxds = TenxDataStorage(sampleid, species=species)
    # tenx.bus_finalize()
    open(finished,"w").write("Completed")

def RunPseudo(sampleid, workflow, full=False):
    workflow.transform (
        name = "kallisto",
        func = RunKallisto,
        args = (
            sampleid,
            pypeliner.managed.OutputFile("kallisto.complete"),
        )
    )
    return workflow
