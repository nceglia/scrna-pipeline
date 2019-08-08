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

# def DownloadFastqs(sampleid, finished):
#     fastqs = [FastQDirectory(fastq, config.prefix, config.jobpath, config.datapath) for fastq in [sampleid]]
#     fastqs = glob.glob(os.path.join(config.datapath,"*.fastq.gz"))
#     assert len(fastqs) > 0, "No Fastqs Download or Found."
#     open(finished,"w").write("Completed")

def RunKallisto(sampleid, finished):
    fastqs = [FastQDirectory(sampleid, config.prefix, config.jobpath, config.datapath)]
    krunner = Kallisto(fastq_directory, sample)
    tenx_path = krunner.count()

def RunBusUpload(sampleid, finished, species):
    print("Uploading ",species, sampleid)
    tenx_output = os.path.join(config.jobpath,".cache/{}/kallisto".format(self.sampleid))
    tenx = TenxAnalysis(tenx_output)
    tenxds = TenxDataStorage(sampleid, species=species)
    print("Running upload")
    tenx.bus_finalize()
    tenxds.upload_kallisto(tenx)
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
    workflow.transform (
        name = "upload_kallisto",
        func = RunBusUpload,
        args = (
            sampleid,
            pypeliner.managed.OutputFile("kallisto.complete"),
        )
    )
    return workflow
