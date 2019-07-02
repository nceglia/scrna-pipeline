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
    for sample in sampleids:
        tenx = TenxDataStorage(sample)
        path = tenx.download()
        path_json = {sample: path}
        open(finished(sample),"w").write(json.dumps(path_json))

def RunExtract(sample_to_path, finished):
    sample = json.loads(open(sample_to_path,"r").read())
    sampleid, path = sample.items().pop()
    tenx_analysis = TenxAnalysis(path)
    tenx_analysis.load()
    tenx_analysis.extract()
    print("Extracted.")
    # qc = QualityControl(tenx_analysis,sampleid)
    # if not os.path.exists(qc.sce):
    #     qc.run(mito=config.mito)



def RunCollection(workflow):
    workflow.transform (
        name = "download_collection",
        func = RunDownload,
        args = (
            open(config.samples, "r").read().splitlines(),
            pypeliner.managed.TempOutputFile("sample_path.json","sample")
        )
    )
    workflow.transform (
        name = "extract_rdata",
        func = RunExtract,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("sample_path.json","sample"),
            pypeliner.managed.OutputFile("sample.rdata","sample")
        )
    )

    return workflow
