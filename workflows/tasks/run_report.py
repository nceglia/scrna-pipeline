import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os

from interface.tenxanalysis import TenxAnalysis
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl
from utils.reporting import Results
from utils.export import exportMD, exportUpload

from utils.config import Configuration, write_config

config = Configuration()

def Run(sampleid, before, finished):
    tenx = TenxDataStorage(sampleid, version="v3")
    tenx.download()
    tenx_analysis = TenxAnalysis(tenx.tenx_path)
    tenx_analysis.load()
    tenx_analysis.extract()
    qc = QualityControl(tenx_analysis,sampleid)
    plots = qc.plots
    cellassign = os.path.join(os.path.split(plots)[0],"cellassignanalysis")
    results = Results(config.jobpath)
    filtered_sce = os.path.join(os.path.split(qc.sce)[0],"sce_cas.rdata")

    results.add_analysis(tenx.tenx_path)
    results.add_sce(filtered_sce)

    #QC Figures
    umi = os.path.join(plots,"umi.png")
    mito = os.path.join(plots,"mito.png")
    ribo = os.path.join(plots, "ribo.png")
    total_counts = os.path.join(plots, "total_counts.png")
    tfbc = os.path.join(plots, "total_features_by_counts.png")
    tcvfc = os.path.join(plots, "total_counts_v_features_by_counts.png")

    #Cellassign
    celltypes = os.path.join(cellassign, "cell_types.png")
    tsne_celltypes = os.path.join(cellassign, "tsne_by_cell_type.png")
    umap_celltypes = os.path.join(cellassign, "umap_by_cell_type.png")

    #Clusters
    tsne_clusters = os.path.join(cellassign, "tsne_by_clusters.png")
    umap_clusters = os.path.join(cellassign, "umap_by_clusters.png")

    results.add_plot(umi,"UMI Distribution")
    results.add_plot(mito,"Mito Distribution")
    results.add_plot(ribo,"Ribo Distribution")
    results.add_plot(total_counts,"Total Counts Distribution")
    results.add_plot(tcvfc,"Total Counts")
    results.add_plot(tcvfc,"Total Features by Counts")

    results.add_plot(tsne_celltypes,"TSNE - Cell Types")
    results.add_plot(umap_celltypes,"UMAP - Cell Types")
    results.add_plot(celltypes,"Cell Types")

    results.add_plot(tsne_clusters,"TSNE - Clusters")
    results.add_plot(umap_clusters,"UMAP - Clusters")

    exportUpload(results)
    open(finished,"w").write("Completed")


def RunReport(sampleid, workflow):
    workflow.transform (
        name = "report",
        func = Run,
        args = (
            sampleid,
            pypeliner.managed.InputFile("cellassignanalysis.complete"),
            pypeliner.managed.OutputFile("report.complete")
        )
    )
    return workflow

if __name__ == '__main__':
    sampleid = sys.argv[1]
    Run(sampleid, None, "qc.complete")
