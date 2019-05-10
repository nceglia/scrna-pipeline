import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys

from interface.tenxanalysis import TenxAnalysis
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl

from utils.config import Configuration, write_config

config = Configuration()

def Run(sampleid, before, finished):
    print("Running QC.")
    tenx = TenxDataStorage(sampleid, version="v3")
    tenx.download()
    tenx_analysis = TenxAnalysis(tenx.tenx_path)
    tenx_analysis.load()
    tenx_analysis.extract()
    print("Extracted.")
    qc = QualityControl(tenx_analysis,sampleid)
    if not os.path.exists(qc.sce):
        qc.run(mito=config.mito)
        print ("Uploading")
        qc.upload_raw()
        qc.upload()
    open(finished,"w").write("Completed")

def RunQC(tenx, workflow):
    workflow.transform (
        name = "quality_control",
        func = Run,
        args = (
            tenx,
            pypeliner.managed.InputFile("upload.complete"),
            pypeliner.managed.OutputFile("qc.complete")
        )
    )
    return workflow
