import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os
import json

from interface.tenxanalysis import TenxAnalysis
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl
from utils.cloud import SampleCollection

from utils.config import Configuration, write_config

config = Configuration()

def RunDownload(sampleids, finished):
    print("Getting Collection.")
    tenx_collection = SampleCollection(sampleids)
    paths = tenx_collection.sample_paths()
    tenx_str = json.loads(str(paths))
    open("sample_paths.json","w").write(json.dumps(tenx_str))
    open(finished,"w").write("Completed")

def RunRdata(sampleid, finished):
    print("Getting Collection.")
    sampleids = open(config.samples, "r").read().splitlines()
    tenx_collection = SampleCollection(sampleids)
    open("sample_paths.json","w").write(json.dumps(tenx_str))
    # for sampleid, tenx_path in tenx_collection:
    #     tenx_analysis = TenxAnalysis(tenx_path)
    #     tenx_analysis.load()
    #     tenx_analysis.extract()
    #     print("Extracted.")
    #     qc = QualityControl(tenx_analysis,sampleid)
    #     if not os.path.exists(qc.sce):
    #         qc.run(mito=config.mito)
    #         print ("Uploading")
    open(finished,"w").write("Completed")


def RunCollection(workflow):
    workflow.transform (
        name = "download_collection",
        func = RunDownload,
        args = (
            open(config.samples, "r").read().splitlines(),
            pypeliner.managed.TempOutputFile("sample_path.json","sample")
        )
    )
    return workflow
