import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os
import json
import shutil
import subprocess
import warnings

from interface.tenxanalysis import TenxAnalysis
from interface.fastqdirectory import FastQDirectory
from utils.cloud import FastqDataStorage
from interface.qualitycontrol import QualityControl
from utils.cloud import SampleCollection
from interface.genemarkermatrix import GeneMarkerMatrix
from utils.plotting import celltypes, tsne_by_cell_type, umap_by_cell_type
from software.cellassign import CellAssign
from software.kallisto import Kallisto

from utils.config import Configuration, write_config

warnings.filterwarnings("ignore")

config = Configuration()

def RunDownload(sampleids, finished):
    for i, sample in enumerate(sampleids):
        fastq = FastqDataStorage(sample)
        fastq.set_data_path(os.path.join(config.datapath,sample))
        path = fastq.download_fastqs()
        path_json = {sample: path}
        open(finished(i),"w").write(json.dumps(path_json))

def RunKallisto(sample_to_path, bus_output):
    sample = json.loads(open(sample_to_path,"r").read())
    sampleid, path = list(sample.items()).pop()
    fastqs = [FastQDirectory(sampleid, config.prefix, config.jobpath, path)]
    krunner = Kallisto(fastqs[0], sampleid)
    bus_path = krunner.count()
    open(bus_output,"w").write(bus_path)


def RunQC(bus_output, sce, filtered_sce):
    bus_path = open(bus_output,"r").read().strip()
    rcode = """
    library(SingleCellExperiment)
    library(scater)
    library(DropletUtils)
    library(scran)
    library(stringr)

    sce <- read10xCounts('{bus_path}')
    sce <- sce[,colSums(counts(sce))>0]
    sce <- sce[rowSums(counts(sce))>0,]
    rowData(sce)$ensembl_gene_id <- rownames(sce)
    sce <- sce[,colSums(counts(sce) != 0) > 100]
    sce <- sce[rowSums(counts(sce))>0,]
    counts(sce) <- data.matrix(counts(sce))

    print("calculating metrics")
    mitochondrial <- as.character(rowData(sce)$Symbol[str_detect(rowData(sce)$Symbol, "^MT\\\-")])
    ribosomal <- as.character(rowData(sce)$Symbol[str_detect(rowData(sce)$Symbol, "^RP(L|S)")])
    rownames(sce) <- rowData(sce)$Symbol
    sce <- calculateQCMetrics(sce, exprs_values = "counts", feature_controls = list(mitochondrial=mitochondrial, ribosomal=ribosomal))
    saveRDS(sce, file='{raw}')
    print("normalizing")
    sce <- normalize(sce)
    print("filtering")
    cells_to_keep <- sce$pct_counts_mitochondrial < 25 && sce$pct_counts_ribosomal < 65
    sce <- sce[,cells_to_keep]
    qclust <- quickCluster(sce, min.size = 100)
    sce <- computeSumFactors(sce, clusters = qclust)
    sce$size_factor <- sizeFactors(sce)
    print("Running red dim")
    sce <- runPCA(sce, ntop = 1000, ncomponents = 50, exprs_values = "logcounts")
    sce <- runTSNE(sce, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)
    sce <- runUMAP(sce, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)
    saveRDS(sce, file='{filtered}')
    """.format(raw=sce,filtered=filtered_sce, bus_path=os.path.abspath(bus_path))
    path = "/".join(bus_path.split("/")[:-1])
    qc_script = os.path.join(path,"convert.R")
    output = open(qc_script,"w")
    output.write(rcode)
    output.close()
    subprocess.call(["Rscript",qc_script])
    output = "/".join(bus_path.split("/")[:-1])
    output = os.path.join(output, "sce.rdata")
    shutil.copyfile(filtered_sce,output)

def RunCellAssign(sce, annot_sce):
    _rho_csv = os.path.join(os.path.split(sce)[0],"rho_csv_sub.csv")
    _fit = os.path.join(os.path.split(sce)[0],"fit_sub.pkl")
    sampleid = sce.split("/")[-2]
    CellAssign.run(sce, config.rho_matrix, _fit, rho_csv=_rho_csv)
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

def RunIntegration(seurats, integrated_seurat):
    rdata = os.path.join(os.path.split(integrated_seurat)[0],"integrate_seurat_cached.rdata")
    object_list = []
    rcode = """
    library(Seurat)
    """
    for idx, object in seurats.items():
        seurat_obj = "seurat{}".format(idx)
        object_list.append(seurat_obj)
        load = """
    seurat{id} <- readRDS("{object}")
        """.format(id=idx,object=object)
        rcode += load
    rcode += """
    object_list <- c({object_list})
    features <- SelectIntegrationFeatures(object.list = c({object_list}), nfeatures = 3000)
    prepped <- PrepSCTIntegration(object.list = c({object_list}), anchor.features = features)
    anchors <- FindIntegrationAnchors(object.list = prepped, normalization.method="SCT", anchor.features=features)
    integrated <- IntegrateData(anchorset = anchors, normalization="SCT")
    saveRDS(integrated, file = "{rdata}")
    integrated <- RunPCA(integrated, verbose = FALSE)
    integrated <- RunUMAP(integrated, dims = 1:30)
    saveRDS(integrated, file ="{rdata}")
    """
    integrate_script = os.path.join(".cache/integration.R")
    output = open(integrate_script,"w")
    output.write(rcode.format(seurat=seurat, object_list=",".join(object_list), rdata=rdata))
    output.close()
    subprocess.call(["Rscript","{}".format(integrate_script)])
    shutil.copyfile(rdata, integrated_seurat)

def RunReport(samples, sces, seurats, tsnes, umaps, tsnecelltypes, umapcelltypes, ridge, features, markers, umi, ribo, mito, counts, raw_sces, summary_path, celltypes, tsne_basic, umap_basic):
    for id, rdata in seurats.items():
        sample_json_path = samples[id]
        sample_json = json.loads(open(sample_json_path,"r").read())
        sample_name = list(sample_json.keys())[0]
        dir = "results"
        if not os.path.exists(dir):
            os.makedirs(dir)

        inventory = dict()
        inventory["seurat"] = os.path.join(dir,"{}_seurat.rds".format(sample_name))
        inventory["sce"] = os.path.join(dir,"{}_sce.rds".format(sample_name))
        inventory["tsne"] = os.path.join(dir,"{}_tsne.png".format(sample_name))
        inventory["umap"] =  os.path.join(dir,"{}_umap.png".format(sample_name))
        inventory["ridge"] = os.path.join(dir,"{}_ridge.png".format(sample_name))
        inventory["features"] = os.path.join(dir,"{}_features.png".format(sample_name))
        inventory["markers"] = os.path.join(dir,"{}_markers.csv".format(sample_name))
        inventory["umi"] =  os.path.join(dir,"{}_umi.png".format(sample_name))
        inventory["ribo"] = os.path.join(dir,"{}_ribo.png".format(sample_name))
        inventory["mito"] = os.path.join(dir,"{}_mito.png".format(sample_name))
        inventory["counts"] = os.path.join(dir,"{}_counts.png".format(sample_name))
        inventory["raw_sce"] = os.path.join(dir,"{}_raw_sce.rdata".format(sample_name))
        inventory["summary"] = os.path.join(dir,"{}_cellranger.html".format(sample_name))
        inventory["tsne_basic"] = os.path.join(dir,"{}_tsne_basic.png".format(sample_name))
        inventory["umap_basic"] = os.path.join(dir,"{}_umap_basic.png".format(sample_name))
        inventory["celltypes"] = os.path.join(dir,"{}_umap_basic.png".format(sample_name))

        shutil.copyfile(rdata, os.path.join(dir,"{}_seurat.rds".format(sample_name)))
        shutil.copyfile(sces[id], os.path.join(dir,"{}_sce.rds".format(sample_name)))
        shutil.copyfile(tsnes[id], os.path.join(dir,"{}_tsne.png".format(sample_name)))
        shutil.copyfile(umaps[id], os.path.join(dir,"{}_umap.png".format(sample_name)))
        shutil.copyfile(tsnecelltypes[id], os.path.join(dir,"{}_tsne_celltype.png".format(sample_name)))
        shutil.copyfile(umapcelltypes[id], os.path.join(dir,"{}_umap_celltype.png".format(sample_name)))
        shutil.copyfile(ridge[id], os.path.join(dir,"{}_ridge.png".format(sample_name)))
        shutil.copyfile(features[id], os.path.join(dir,"{}_features.png".format(sample_name)))
        shutil.copyfile(markers[id], os.path.join(dir,"{}_markers.csv".format(sample_name)))
        shutil.copyfile(umi, os.path.join(dir,"{}_umi.png".format(sample_name)))
        shutil.copyfile(ribo, os.path.join(dir,"{}_ribo.png".format(sample_name)))
        shutil.copyfile(mito, os.path.join(dir,"{}_mito.png".format(sample_name)))
        shutil.copyfile(counts, os.path.join(dir,"{}_counts.png".format(sample_name)))
        shutil.copyfile(raw_sces, os.path.join(dir,"{}_raw_sce.rdata".format(sample_name)))
        shutil.copyfile(summary_path[id], os.path.join(dir,"{}_cellranger.html".format(sample_name)))
        shutil.copyfile(tsne_basic, inventory["tsne_basic"])
        shutil.copyfile(tsne_basic, inventory["umap_basic"])
        shutil.copyfile(tsne_basic, inventory["celltypes"])

        generateHTML(dir, sample_name, inventory)

        # report = ReportStorage(dir)
        # report.upload(os.path.dirname(os.path.realpath(__file__)), sample_name)



def RunPipeline(workflow):
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
        name = "run_kallisto",
        func = RunKallisto,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("sample_path.json","sample"),
            pypeliner.managed.TempOutputFile("bus_path.json","sample")
        )
    )

    workflow.transform (
        name = "run_qc",
        func = RunQC,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("bus_path.json","sample"),
            pypeliner.managed.TempOutputFile("sce.rdata","sample"),
            pypeliner.managed.TempOutputFile("filtered.rdata","sample"),
        )
    )

    workflow.transform (
        name = "run_cellassign",
        func = RunCellAssign,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("filtered.rdata","sample"),
            pypeliner.managed.TempOutputFile("annot.rdata","sample"),
        )
    )

    workflow.transform (
        name = "run_convert",
        func = RunConvert,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("annot.rdata","sample"),
            pypeliner.managed.TempOutputFile("seurat.rdata","sample"),
        )
    )

    workflow.transform (
        name = "run_seurat",
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
        name = "integrate",
        func = RunIntegration,
        args = (
            pypeliner.managed.TempInputFile("seurat_qcd.rdata","sample"),
            pypeliner.managed.TempOutputFile("seurat_integrated.rdata")
        )
    )

    # workflow.transform (
    #     name = "report",
    #     func = RunReport,
    #     args = (
    #         pypeliner.managed.TempInputFile("sample_path.json","sample"),
    #         pypeliner.managed.TempInputFile("sce.rdata","sample"),
    #         pypeliner.managed.TempInputFile("seurat_qcd.rdata","sample"),
    #         pypeliner.managed.TempInputFile("seurat_tsne.png","sample"),
    #         pypeliner.managed.TempInputFile("seurat_umap.png","sample"),
    #         pypeliner.managed.TempInputFile("seurat_tsne_celltype.png","sample"),
    #         pypeliner.managed.TempInputFile("seurat_umap_celltype.png","sample"),
    #         pypeliner.managed.TempInputFile("seurat_ridge.png","sample"),
    #         pypeliner.managed.TempInputFile("seurat_features.png","sample"),
    #         pypeliner.managed.TempInputFile("markers.csv","sample"),
    #         pypeliner.managed.TempInputFile("seurat_integrated.rdata"),
    #         # pypeliner.managed.TempInputFile("umi.png"),
    #         # pypeliner.managed.TempInputFile("mito.png"),
    #         # pypeliner.managed.TempInputFile("ribo.png"),
    #         # pypeliner.managed.TempInputFile("counts.png"),
    #         # pypeliner.managed.TempInputFile("raw_sce.rdata"),
    #         pypeliner.managed.TempInputFile("summary_path.html","sample"),
    #         # pypeliner.managed.TempInputFile("celltypes.png"),
    #         pypeliner.managed.TempInputFile("tsne_by_celltype.png"),
    #         pypeliner.managed.TempInputFile("umap_by_celltype.png"),
    #     )
    # )


    return workflow


# def generateHTML(report, sampleid, inventory):
#     html = os.path.join(report, "report.html")
#     output = open(html,"w")
#     output.write(template_header)
#     res = template.format(sampleid = sampleid,
#                             raw_sce = inventory["raw_sce"],
#                             sce = inventory["sce"],
#                             seurat = inventory["seurat"],
#                             markers = inventory["markers"],
#                             summary = inventory["summary"],
#                             umi = inventory["umi"],
#                             mito = inventory["mito"],
#                             ribo = inventory["ribo"],
#                             counts = inventory["counts"],
#                             tsne = inventory["tsne"],
#                             umap = inventory["umap"],
#                             tsne_celltype_seurat = inventory["tsne_celltype_seurat"],
#                             umap_celltype_seurat = inventory["umap_celltype_seurat"],
#                             tsne_basic = inventory["tsne_basic"],
#                             umap_basic = inventory["umap_basic"],
#                             celltypes = inventory["celltypes"],
#                             ridge = inventory["ridge"],
#                             features = inventory["features"]
#                             )
#     output.write(res)
#     output.write(template_footer)
#     output.close()


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
            <tr><td>Raw SCE:</td><td>{raw_sce}</td></tr>
            <tr><td>SCE (cellassign):</td><td>{sce}</td></tr>
            <tr><td>Seurat (cellassign):</td><td>{seurat}</td></tr>
            <tr><td>Markers by Cluster:</td><td>{markers}</td></tr>
            <tr><td>Cellranger Summary:</td><td>{summary}</td></tr>
        </table>
        <br>
        <h3>QC</h3>
        <table width="80%">
            <tr><td>UMI</td><td>Mito</td><td>Ribo</td><td>Counts</td></tr>
            <tr><td><img src="{umi}"></td><td><img src="{mito}"></td><td><img src="{ribo}"></td><td><img src="{counts}"></td></tr>
        </table>
        <br>
        <h3>Analysis</h3>
        <table width="80%">
          <tr><td>Cell Types</td></tr>
          <tr><td><img src="{celltypes}"></td></tr>

          <tr><td>UMAP</td><td>TSNE</td></tr>
          <tr><td><img src="{umap_basic}"></td><td><img src="{tsne_basic}"></td></tr>

          <tr><td>UMAP CellType (Seurat)</td><td>TSNE Celltype (Seurat)</td></tr>
          <tr><td><img src="{seurat_umap_celltype}"></td><td><img src="{seurat_tsne_celltype}"></td></tr>

          <tr><td>Markers in Clusters</td><td>Markers Embedding</td></tr>
          <tr><td><img src="{ridge}"></td><td><img src="{features}"></td></tr>
        </table>
        <br><br><br>
        <script type="text/javascript"charset="utf-8">
            d3.text("{markers}", function(data) {
                var parsedCSV = d3.csv.parseRows(data);
"""
template_footer = """
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
