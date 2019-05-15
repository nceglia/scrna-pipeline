import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import scanorama
import sys
import os
import pickle
import glob

from interface.tenxanalysis import TenxAnalysis
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl
from software.batchcorrection import Scanorama

from utils.config import Configuration

config = Configuration()

def Search(sampleid):
    tenxs = []
    tenx = TenxDataStorage(sampleid, version="v3")
    tenx.download()
    analysis_path = tenx.tenx_path
    print(analysis_path)
    tenx_analysis = TenxAnalysis(analysis_path)
    tenx_analysis.load()
    tenx_analysis.extract()
    qc = QualityControl(tenx_analysis, sampleid)
    tenxs.append(tenx_analysis.adata(qc.sce))
    print ("Loading main sce {}".format(sampleid))
    sys.stdout.flush()
    samples = glob.glob("../../*/runs/.cache/*/metrics_summary.csv")
    for sample in samples[:2]:
        print ("Loading project sample {}".format(sample))
        sys.stdout.flush()
        sample_rel_path = os.path.split(sample)[0]
        sid = sample_rel_path.split("/")[-1]
        sidsce = os.path.join(sample_rel_path,"{0}.rdata".format(sid))
        if not os.path.exists(sidsce):
            print("Not found",sidsce)
            continue
        tenx_analysis = TenxAnalysis(sample_rel_path)
        tenx_analysis.load()
        tenx_analysis.extract()
        tenxs.append(tenx_analysis.adata(sidsce))
    print ("Finished project tree search.")
    sys.stdout.flush()
    return tenxs

def Run(sampleid, before, finished):
    adatas = Search(sampleid)
    print ("Correcting on {} samples.".format(len(adatas)))
    sys.stdout.flush()
    corrected = Scanorama.correct(adatas)
    sys.stdout.flush()
    if not os.path.exists(".cache/corrected"):
        os.makedirs(".cache/corrected")
    TenxAnalysis.make_10x_output(str(corrected[0]),".cache/corrected")
    open(finished,"w").write("Completed")


def RunCorrection(sampleid, workflow):
    workflow.transform (
        name = "batch_correct",
        func = Run,
        args = (
            sampleid,
            pypeliner.managed.InputFile("qc.complete"),
            pypeliner.managed.OutputFile("correction.complete")
        )
    )
    return workflow

if __name__ == '__main__':
    Run("AS_CD45P_IGO_09443_F_16", "qc.complete", "correction.complete")
