import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os
import json
import shutil
import subprocess

from interface.tenxanalysis import TenxAnalysis
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl
from utils.cloud import SampleCollection
from interface.qualitycontrol import QualityControl
from interface.genemarkermatrix import GeneMarkerMatrix
from utils.plotting import celltypes, tsne_by_cell_type, umap_by_cell_type
from software.cellassign import CellAssign

from utils.config import Configuration, write_config

config = Configuration()

def RunDownload(sampleids, finished):
    for i, sample in enumerate(sampleids):
        tenx = TenxDataStorage(sample)
        path = tenx.download()
        path_json = {sample: path}
        open(finished(i),"w").write(json.dumps(path_json))

def RunExtract(sample_to_path, rdata_path, summmary_path):
    sample = json.loads(open(sample_to_path,"r").read())
    sampleid, path = list(sample.items()).pop()
    tenx_analysis = TenxAnalysis(path)
    tenx_analysis.load()
    tenx_analysis.extract()
    qc = QualityControl(tenx_analysis, sampleid)
    if not os.path.exists(qc.sce):
        qc.run(mito=config.mito)
    shutil.copyfile(tenx_analysis.summary(), summary_path)
    shutil.copyfile(qc.sce, rdata_path)

def RunCellAssign(sce, annot_sce, rho_csv, fit):
    sampleid = sce.split("/")[-2]
    CellAssign.run(sce, config.rho_matrix, fit, rho_csv=rho_csv)
    filtered_sce = os.path.join(os.path.split(sce)[0],"sce_cas.rdata")
    shutil.copyfile(filtered_sce, annot_sce)

def RunConvert(sce, seurat):
    seurat_cached = os.path.join(os.path.split(sce)[0],"seurat_raw.rdata")
    sce_cached = os.path.join(os.path.split(sce)[0],"sce_cas.rdata")
    rcode = """
    library(Seurat)
    library(SingleCellExperiment)
    library(scater)
    sce <- readRDS('{sce}')
    rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ensembl_gene_id, rownames(sce))
    seurat <- as.Seurat(sce, counts = "counts", data = "logcounts")
    saveRDS(seurat,file='{seurat}')
    """
    path = os.path.split(sce)[0]
    convert_script = os.path.join(path,"convert.R")
    output = open(convert_script,"w")
    output.write(rcode.format(sce=sce_cached,seurat=seurat_cached))
    output.close()
    subprocess.call(["Rscript","{}".format(convert_script)])
    shutil.copyfile(seurat_cached, seurat)

def RunSeuratWorkflow(seurat, qcd_seurat):
    seurat_cached = os.path.join(os.path.split(seurat)[0],"seuret_annot.rdata")
    rcode = """
    library(Seurat)
    seurat <- readRDS("{seurat}")
    seurat <- NormalizeData(object = seurat)
    seurat <- FindVariableFeatures(object = seurat)
    seurat <- ScaleData(object = seurat)
    seurat <- RunPCA(object = seurat)
    seurat <- FindNeighbors(object = seurat)
    seurat <- FindClusters(object = seurat)
    seurat <- RunTSNE(object = seurat)
    seurat <- RunUMAP(object = seurat, reduction = "pca", dims = 1:20)
    saveRDS(seurat, file = '{qcd_seurat}')"""
    path = os.path.split(seurat)[0]
    qc_script = os.path.join(path,"qc.R")
    output = open(qc_script,"w")
    output.write(rcode.format(seurat=seurat, qcd_seurat=seurat_cached))
    output.close()
    subprocess.call(["Rscript", "{}".format(qc_script)])
    shutil.copyfile(seurat_cached, qcd_seurat)

def RunSeuratViz(seurat, tsne, umap, tsne_celltype, umap_celltype, ridge, exprs):
    marker_list = GeneMarkerMatrix.read_yaml(config.rho_matrix)
    markers = ["'" + marker + "'" for marker in marker_list.genes]
    tsne_plot = os.path.join(os.path.split(seurat)[0],"tsne.png")
    umap_plot = os.path.join(os.path.split(seurat)[0],"umap.png")
    tsne_celltype_plot = os.path.join(os.path.split(seurat)[0],"tsne_celltype.png")
    umap_celltype_plot = os.path.join(os.path.split(seurat)[0],"umap_celltype.png")
    ridge_plot = os.path.join(os.path.split(seurat)[0],"ridge.png")
    exprs_plot = os.path.join(os.path.split(seurat)[0],"features.png")
    rcode = """
    library(Seurat)
    library(ggplot2)
    seurat <- readRDS("{seurat}")

    png("{tsne}")
    DimPlot(object = seurat, reduction = "tsne")
    dev.off()
    png("{umap}")
    DimPlot(object = seurat, reduction = "umap")
    dev.off()

    png("{tsne_celltype}")
    DimPlot(object = seurat, reduction = "tsne", group.by = "cell_type")
    dev.off()
    png("{umap_celltype}")
    DimPlot(object = seurat, reduction = "umap", group.by = "cell_type")
    dev.off()

    png("{ridge}",width=600,heigh=5000)
    RidgePlot(object = seurat, features = c({markers}), ncol = 2)
    dev.off()

    png("{exprs}",width=600,heigh=5000)
    FeaturePlot(object = seurat, features = c({markers}), ncol= 2)
    dev.off()

    """
    path = os.path.split(seurat)[0]
    qc_script = os.path.join(path,"viz.R")
    output = open(qc_script,"w")
    output.write(rcode.format(seurat=seurat, tsne=tsne_plot, umap=umap_plot, tsne_celltype=tsne_celltype_plot, umap_celltype=umap_celltype_plot, markers=",".join(markers), ridge = ridge_plot, exprs=exprs_plot))
    output.close()
    subprocess.call(["Rscript","{}".format(qc_script)])
    shutil.copyfile(tsne_plot, tsne)
    shutil.copyfile(umap_plot, umap)
    shutil.copyfile(tsne_celltype_plot, tsne_celltype)
    shutil.copyfile(umap_celltype_plot, umap_celltype)
    shutil.copyfile(ridge_plot, ridge)
    shutil.copyfile(exprs_plot, exprs)

def RunMarkers(seurat,marker_table):
    marker_csv_cached = os.path.join(os.path.split(seurat)[0],"marker_table.csv")
    rcode = """
    library(Seurat)
    library(dplyr)
    seurat <- readRDS("{seurat}")
    markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    marker_table <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
    marker_table <- as.data.frame(marker_table)
    write.csv(marker_table, file = "{marker_csv}")
    """
    path = os.path.split(seurat)[0]
    marker_script = os.path.join(path,"markers.R")
    output = open(marker_script,"w")
    output.write(rcode.format(seurat=seurat, marker_csv=marker_csv_cached))
    output.close()
    subprocess.call(["Rscript","{}".format(marker_script)])
    shutil.copyfile(marker_csv_cached, marker_table)

def RunReport(samples, sces, seurats, tsnes, umaps, tsnecelltypes, umapcelltypes, ridge, exprs, markers, umi, ribo, mito, counts, raw_sces, summary_path):

    for id, rdata in seurats.items():
        sample_json_path = samples[id]
        sample_json = json.loads(open(sample_json_path,"r").read())
        sample_name = list(sample_json.keys())[0]
        dir = "results"
        if not os.path.exists(dir):
            os.makedirs(dir)
        shutil.copyfile(rdata, os.path.join(dir,"{}_seurat.rds".format(sample_name)))
        shutil.copyfile(sces[id], os.path.join(dir,"{}_sce.rds".format(sample_name)))
        shutil.copyfile(tsnes[id], os.path.join(dir,"{}_tsne.png".format(sample_name)))
        shutil.copyfile(umaps[id], os.path.join(dir,"{}_umap.png".format(sample_name)))
        shutil.copyfile(tsnecelltypes[id], os.path.join(dir,"{}_tsne_celltype.png".format(sample_name)))
        shutil.copyfile(umapcelltypes[id], os.path.join(dir,"{}_umap_celltype.png".format(sample_name)))
        shutil.copyfile(tsnecelltypes[id], os.path.join(dir,"{}_ridge.png".format(sample_name)))
        shutil.copyfile(umapcelltypes[id], os.path.join(dir,"{}_features.png".format(sample_name)))
        shutil.copyfile(markers[id], os.path.join(dir,"{}_markers.csv".format(sample_name)))
        print(umi[id])
        shutil.copyfile(umi[id], os.path.join(dir,"{}_umi.png".format(sample_name)))
        shutil.copyfile(ribo[id], os.path.join(dir,"{}_ribo.png".format(sample_name)))
        shutil.copyfile(mito[id], os.path.join(dir,"{}_mito.png".format(sample_name)))
        shutil.copyfile(raw_sces[id], os.path.join(dir,"{}_raw_sce.rdata".format(sample_name)))
        shutil.copyfile(summary_path[id], os.path.join(dir,"{}_cellranger.html".format(sample_name)))

        report = ReportStorage(dir)
        report.upload(os.path.dirname(os.path.realpath(__file__)), sample_name)



def RunCollection(workflow):
    all_samples = open(config.samples, "r").read().splitlines()
    workflow.transform (
        name = "download_collection",
        func = RunDownload,
        args = (
            all_samples,
            pypeliner.managed.TempOutputFile("sample_path.json","sample")
        )
    )
    workflow.transform (
        name = "extract_rdata",
        func = RunExtract,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("sample_path.json","sample"),
            pypeliner.managed.TempOutputFile("sample.rdata","sample"),
            pypeliner.managed.TempOutputFile("summary_path.json","sample")
        )
    )

    workflow.transform (
        name = "run_cellassign",
        func = RunCellAssign,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("sample.rdata","sample"),
            pypeliner.managed.TempOutputFile("sce.rdata","sample"),
            pypeliner.managed.TempOutputFile("rho.csv","sample"),
            pypeliner.managed.TempOutputFile("cellassignfit.rdata","sample"),
        )
    )

    workflow.transform (
        name = "run_convert",
        func = RunConvert,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("sce.rdata","sample"),
            pypeliner.managed.TempOutputFile("seurat.rdata","sample"),
        )
    )

    workflow.transform (
        name = "run_qc",
        func = RunSeuratWorkflow,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("seurat.rdata","sample"),
            pypeliner.managed.TempOutputFile("seurat_qcd.rdata","sample"),
        )
    )

    workflow.transform (
        name = "visualize_sample",
        func = RunSeuratViz,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("seurat_qcd.rdata","sample"),
            pypeliner.managed.TempOutputFile("seurat_tsne.png","sample"),
            pypeliner.managed.TempOutputFile("seurat_umap.png","sample"),
            pypeliner.managed.TempOutputFile("seurat_tsne_celltype.png","sample"),
            pypeliner.managed.TempOutputFile("seurat_umap_celltype.png","sample"),
            pypeliner.managed.TempOutputFile("seurat_ridge.png","sample"),
            pypeliner.managed.TempOutputFile("seurat_features.png","sample"),
        )
    )

    workflow.transform (
        name = "find_markers",
        func = RunMarkers,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("seurat_qcd.rdata","sample"),
            pypeliner.managed.TempOutputFile("markers.csv","sample"),
        )
    )

    workflow.transform (
        name = "report",
        func = RunReport,
        args = (
            pypeliner.managed.TempInputFile("sample_path.json","sample"),
            pypeliner.managed.TempInputFile("sce.rdata","sample"),
            pypeliner.managed.TempInputFile("seurat_qcd.rdata","sample"),
            pypeliner.managed.TempInputFile("seurat_tsne.png","sample"),
            pypeliner.managed.TempInputFile("seurat_umap.png","sample"),
            pypeliner.managed.TempInputFile("seurat_tsne_celltype.png","sample"),
            pypeliner.managed.TempInputFile("seurat_umap_celltype.png","sample"),
            pypeliner.managed.TempInputFile("seurat_ridge.png","sample"),
            pypeliner.managed.TempInputFile("seurat_features.png","sample"),
            pypeliner.managed.TempInputFile("markers.csv","sample"),
            pypeliner.managed.TempInputFile("umi.png"),
            pypeliner.managed.TempInputFile("mito.png"),
            pypeliner.managed.TempInputFile("ribo.png"),
            pypeliner.managed.TempInputFile("counts.png"),
            pypeliner.managed.TempInputFile("raw_sce.rdata"),
            pypeliner.managed.TempInputFile("summary_path.html","sample")
        )
    )


    return workflow



template_header = """
<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="utf-8">
        <style>
            table {
                border-collapse: collapse;
                border: 2px black solid;
                font: 12px sans-serif;
            }

            td {
                border: 1px black solid;
                padding: 5px;
            }
        </style>
    </head>
    <body>
        <script src="http://d3js.org/d3.v3.min.js"></script>
        <!-- <script src="d3.min.js?v=3.2.8"></script> -->
        <center>"""
template = """
        <h1><font face="helvetica" color="#1F1F1F">{sampleid}</font></h1><br><br>
        <h3>Summary</h3>
        <br>
        <table>
            <tr><td>Cells:</td><td>{cells}</td></tr>
            <tr><td>Raw SCE:</td><td>{raw_sce}</td></tr>
            <tr><td>SCE (cellassign):</td><td>{sce}</td></tr>
            <tr><td>Seurat (cellassign):</td><td>{seurat}</td></tr>
            <tr><td>Markers by Cluster:</td><td>{markers}</td></tr>
            <tr><td>Cellranger Summary:</td><td>{summary}</td></tr>
        </table>
        <br>
        <table>
            <tr><td>Cells:</td><td>{cells}</td></tr>
            <tr><td>Raw SCE:</td><td>{raw_sce}</td></tr>
            <tr><td>SCE (cellassign):</td><td>{sce}</td></tr>
            <tr><td>Seurat (cellassign):</td><td>{seurat}</td></tr>
            <tr><td>Cellranger Summary:</td><td>{summary}</td></tr>
        </table>
        <br>
        <table width="80%">
          <tr><td>UMAP CellType</td><td>TSNE Celltype</td></tr>
          <tr><td><img src="umap_celltype.png"></td><td><img src="tsne_celltype.png"></td></tr>
          <tr><td>UMAP</td><td>TSNE</td></tr>
          <tr><td><img src="umap.png"></td><td><img src="tsne.png"></td></tr>
          <tr><td>Markers in Clusters</td><td>Markers Embedding</td></tr>
          <tr><td><img src="ridge_markers.png"></td><td><img src="feature_markers.png"></td></tr>
        </table>
"""
template_footer = """
        <br><br><br>
        <script type="text/javascript"charset="utf-8">
            d3.text("marker_table.csv", function(data) {
                var parsedCSV = d3.csv.parseRows(data);

                var container = d3.select("body")
                    .append("table").attr("width","100%")

                    .selectAll("tr")
                        .data(parsedCSV).enter()
                        .append("tr")

                    .selectAll("td")
                        .data(function(d) { return d; }).enter()
                        .append("td")
                        .text(function(d) { return d; });
            });
        </script>
        </center>
    </body>
</html>
"""
