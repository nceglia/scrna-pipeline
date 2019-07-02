import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os
import json
import shutil

from interface.tenxanalysis import TenxAnalysis
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl
from utils.cloud import SampleCollection

from utils.config import Configuration, write_config

config = Configuration()

def RunDownload(sampleids, finished):
    print(sampleids)
    for sample in sampleids:
        tenx = TenxDataStorage(sample)
        path = tenx.download()
        path_json = {sample: path}
        open(finished(sample),"w").write(json.dumps(path_json))

def RunExtract(sample_to_path, rdata_path):
    sample = json.loads(open(sample_to_path,"r").read())
    print(sample)
    sampleid, path = list(sample.items()).pop()
    tenx_analysis = TenxAnalysis(path)
    tenx_analysis.load()
    tenx_analysis.extract()
    qc = QualityControl(tenx_analysis,sampleid)
    if not os.path.exists(qc.sce):
        qc.run(mito=config.mito)
    shutil.copyfile(qc.sce, rdata_path)


def RunCollect(rdata, manifest):
    output = open(manifest,"w")
    for id, rdata in rdata.items():
        output.write(rdata+"\n")
    output.close()

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

    workflow.transform (
        name = "create_manifest",
        func = RunCollect,
        axes = ('sample',),
        args = (
            pypeliner.managed.InputFile("sample.json","sample"),
            pypeliner.managed.OutputFile("manifest.txt")
        )
    )

    return workflow