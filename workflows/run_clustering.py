import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os
import pickle

from interface.tenxanalysis import TenxAnalysis
from software.cellassign import CellAssign
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl
from utils.plotting import tsne_by_cluster, umap_by_cluster


from utils.config import Configuration

def Run(sampleid, before, finished):
    clustering = ".cache/{}/clustering/".format(sampleid)
    if not os.path.exists(clustering):
        os.makedirs(clustering)
    cluster_results = os.path.join(clustering, "{}_clusters.pkl".format(sampleid))
    if not os.path.exists(cluster_results):
        tenx = TenxDataStorage(sampleid, version="v3")
        tenx.download()
        analysis_path = tenx.tenx_path
        tenx_analysis = TenxAnalysis(analysis_path)
        tenx_analysis.load()
        tenx_analysis.extract()
        qc = QualityControl(tenx_analysis, sampleid)
        clusters = tenx_analysis.clusters(qc.sce)
        pickle.dump(clusters, open(cluster_results,"wb"))
    else:
        clusters = pickle.load(open(cluster_results,"rb"))
    tsne_by_cluster(qc.sce, clusters, sampleid, clustering)
    umap_by_cluster(qc.sce, clusters, sampleid, clustering)
    open(finished,"w").write("Completed")


def RunClustering(sampleid, workflow):
    workflow.transform (
        name = "clustering",
        func = Run,
        args = (
            sampleid,
            pypeliner.managed.InputFile("cellassignanalysis.complete"),
            pypeliner.managed.OutputFile("clustering.complete")
        )
    )
    return workflow
