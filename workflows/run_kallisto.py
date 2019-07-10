import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os

from software.kallisto import Kallisto
from interface.fastqdirectory import FastQDirectory
from interface.tenxanalysis import TenxAnalysis
from utils.cloud import TenxDataStorage
from utils.config import Configuration, write_config

config = Configuration()

def RunKallisto(sampleid, finished):
    fastqs = [FastQDirectory(sampleid, config.prefix, config.jobpath, config.datapath)]
    krunner = Kallisto(fastq_directory, tenx_analysis)

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
