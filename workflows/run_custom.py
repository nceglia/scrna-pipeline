import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os
import json
import shutil
import subprocess
import collections
import numpy
import pickle
import pyparsing as pp

from interface.tenxanalysis import TenxAnalysis
from interface.singlecellexperiment import SingleCellExperiment
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl
from utils.cloud import SampleCollection, VizReportStorage
from interface.qualitycontrol import QualityControl
from interface.genemarkermatrix import GeneMarkerMatrix
from utils.plotting import celltypes, tsne_by_cell_type, umap_by_cell_type
from software.cellassign import CellAssign

from utils.config import Configuration, write_config

config = Configuration()

def RunParse(sampleids, finished):
    if not os.path.exists(os.path.join(config.jobpath,"results")):
        os.makedirs(os.path.join(config.jobpath,"results"))
    for i, sample in enumerate(sampleids):
        path = os.path.join(config.datapath,sample)
        path_json = {sample: path}
        open(finished(i),"w").write(json.dumps(path_json))

def RunQC(custom_output, sce, filtered_sce):
    sample = json.loads(open(custom_output,"r").read())
    sampleid, path = list(sample.items()).pop()
    cached_sce = os.path.join(config.jobpath,"results","sce_{}.rdata".format(sampleid))
    rcode = """
    library(SingleCellExperiment)
    library(scater)
    library(DropletUtils)
    library(scran)
    library(stringr)

    sce <- read10xCounts('{path}')
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
    cells_to_keep <- sce$pct_counts_mitochondrial < 25
    sce <- sce[,cells_to_keep]
    cells_to_keep <- sce$pct_counts_ribosomal < 65
    sce <- sce[,cells_to_keep]
    qclust <- quickCluster(sce, min.size = 100)
    sce <- computeSumFactors(sce, clusters = qclust)
    sce$size_factor <- sizeFactors(sce)
    print("Running red dim")
    sce <- runPCA(sce, ntop = 1000, ncomponents = 50, exprs_values = "logcounts")
    sce <- runTSNE(sce, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)
    sce <- runUMAP(sce, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)
    saveRDS(sce, file='{filtered}')
    """.format(raw=sce,filtered=cached_sce, path=os.path.abspath(path))
    path = os.path.split(sce)[0]
    qc_script = os.path.join(path,"convert_{}.R".format(sampleid))
    output = open(qc_script,"w")
    output.write(rcode)
    output.close()
    if not os.path.exists(cached_sce):
        subprocess.call(["Rscript",qc_script])
    shutil.copyfile(cached_sce,filtered_sce)
    shutil.copyfile(cached_sce,sce)

def RunCellAssign(custom_output, sce, annot_sce):
    sample = json.loads(open(custom_output,"r").read())
    sampleid, path = list(sample.items()).pop()
    temp = os.path.split(annot_sce)[0]
    _rho_csv = os.path.join(os.path.split(sce)[0],"rho_csv_sub.csv")
    _fit = os.path.join(os.path.split(sce)[0],"fit_sub.pkl")
    _filtered_sce = os.path.join(os.path.split(sce)[0],"sce_cas.rdata")
    filtered_sce = os.path.join(config.jobpath,"results","sce_cas_{}.rdata".format(sampleid))
    if not os.path.exists(filtered_sce):
        rho = os.path.join(config.jobpath,config.rho_matrix)
        CellAssign.run(sce, rho, _fit, rho_csv=_rho_csv,lsf=True)
        shutil.copyfile(_filtered_sce, filtered_sce)
    shutil.copyfile(filtered_sce, annot_sce)

def RunConvert(custom_output, sce, seurat):
    sample = json.loads(open(custom_output,"r").read())
    sampleid, path = list(sample.items()).pop()
    seurat_cached = os.path.join(config.jobpath,"results","{}_seurat_raw.rdata".format(sampleid))
    # sce_cached = os.path.join(config.jobpath,"results","sce_cas_{}.rdata".format(sampleid))
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
    convert_script = os.path.join(path,"convert_{}.R".format(sampleid))
    output = open(convert_script,"w")
    output.write(rcode.format(sce=sce,seurat=seurat_cached))
    output.close()
    if not os.path.exists(seurat_cached):
        subprocess.call(["Rscript","{}".format(convert_script)])
    shutil.copyfile(seurat_cached, seurat)

def RunSeuratWorkflow(custom_output, seurat, qcd_seurat, qcd_sce):
    sample = json.loads(open(custom_output,"r").read())
    sampleid, path = list(sample.items()).pop()
    seurat_cached = os.path.join(config.jobpath, "results","seuret_annot_{}.rdata".format(sampleid))
    sce_cached = os.path.join(config.jobpath, "results","sce_annot_{}.rdata".format(sampleid))
    umap = os.path.join(config.jobpath, "results","umap_cluster_{}.png".format(sampleid))
    tsne = os.path.join(config.jobpath, "results","tsne_cluster_{}.png".format(sampleid))
    rcode = """
    library(Seurat)
    library(sctransform)
    library(SingleCellExperiment)
    seurat <- readRDS("{seurat}")
    seurat <- SCTransform(object = seurat)
    seurat <- RunPCA(object = seurat)
    seurat <- FindNeighbors(object = seurat)
    seurat <- FindClusters(object = seurat)
    seurat <- RunTSNE(object = seurat, reduction = "pca")
    seurat <- RunUMAP(object = seurat, reduction = "pca", dims = 1:50)
    png("{tsne}")
    DimPlot(object = seurat, reduction = "tsne")
    dev.off()
    png("{umap}")
    DimPlot(object = seurat, reduction = "umap")
    dev.off()
    saveRDS(seurat, file = '{qcd_seurat}')
    sce <- as.SingleCellExperiment(seurat)
    rowData(sce)$Symbol <- rownames(sce)
    saveRDS(sce, file="{qcd_sce}")
    """
    path = os.path.split(seurat)[0]
    qc_script = os.path.join(path,"qc.R")
    output = open(qc_script,"w")
    output.write(rcode.format(seurat=seurat, qcd_seurat=seurat_cached, qcd_sce=sce_cached, umap=umap, tsne=tsne))
    output.close()
    if not os.path.exists(seurat_cached) or not os.path.exists(sce_cached):
        subprocess.call(["Rscript", "{}".format(qc_script)])
    shutil.copyfile(seurat_cached, qcd_seurat)
    shutil.copyfile(sce_cached, qcd_sce)

def RunSeuratViz(custom_output, seurat, tsne_celltype, umap_celltype, ridge, exprs):
    sample = json.loads(open(custom_output,"r").read())
    sampleid, path = list(sample.items()).pop()
    marker_list = GeneMarkerMatrix.read_yaml(config.rho_matrix)
    markers = ["'" + marker + "'" for marker in marker_list.genes]
    tsne_celltype_plot = os.path.join(config.jobpath,"results","tsne_celltype_{}_sctransform.png".format(sampleid))
    umap_celltype_plot = os.path.join(config.jobpath,"results","umap_celltype_{}_sctransform.png".format(sampleid))
    ridge_plot = os.path.join(config.jobpath,"results","ridge_{}.png".format(sampleid))
    exprs_plot = os.path.join(config.jobpath,"results","features_{}.png".format(sampleid))
    rcode = """
    library(Seurat)
    library(ggplot2)
    seurat <- readRDS("{seurat}")

    png("{umap_celltype}")
    DimPlot(object = seurat, reduction = "umap", group.by = "cell_type")
    dev.off()

    png("{tsne_celltype}")
    DimPlot(object = seurat, reduction = "tsne", group.by = "cell_type")
    dev.off()

    png("{ridge}",width=600,height=5000)
    RidgePlot(object = seurat, features = c({markers}), ncol = 2)
    dev.off()

    png("{exprs}",width=600,height=5000)
    FeaturePlot(object = seurat, features = c({markers}), ncol= 2)
    dev.off()
    """
    path = os.path.split(seurat)[0]
    qc_script = os.path.join(path,"viz_{}.R".format(sampleid))
    output = open(qc_script,"w")
    output.write(rcode.format(seurat=seurat, tsne_celltype=tsne_celltype_plot, umap_celltype=umap_celltype_plot, markers=",".join(markers), ridge = ridge_plot, exprs=exprs_plot))
    output.close()
    if not os.path.exists(exprs_plot):
        subprocess.call(["Rscript","{}".format(qc_script)])
    shutil.copyfile(tsne_celltype_plot, tsne_celltype)
    shutil.copyfile(umap_celltype_plot, umap_celltype)
    shutil.copyfile(ridge_plot, ridge)
    shutil.copyfile(exprs_plot, exprs)

def RunMarkers(custom_output,seurat,marker_table):
    sample = json.loads(open(custom_output,"r").read())
    sampleid, path = list(sample.items()).pop()
    marker_csv_cached = os.path.join(config.jobpath, "results","marker_table_{}.csv".format(sampleid))
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
    marker_script = os.path.join(path,"markers_{}.R".format(sampleid))
    output = open(marker_script,"w")
    output.write(rcode.format(seurat=seurat, marker_csv=marker_csv_cached))
    output.close()
    if not os.path.exists(marker_csv_cached):
        subprocess.call(["Rscript","{}".format(marker_script)])
    shutil.copyfile(marker_csv_cached, marker_table)

def RunSeuratIntegration(seurats, integrated_seurat, integrated_sce, integrated_tsne, integrated_umap):
    rdata = os.path.join(config.jobpath,"results","integrated_seurat_seurat.rdata")
    sce_cached = os.path.join(config.jobpath,"results","integrated_seurat_sce.rdata")
    umap = os.path.join(config.jobpath,"results","integrated_seurat_tsne.rdata")
    tsne = os.path.join(config.jobpath,"results","integrated_seurat_umap.rdata")
    object_list = []
    rcode = """
    library(Seurat)
    library(SingleCellExperiment)
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
    sce <- as.SingleCellExperiment(integrated)
    rowData(sce)$Symbol <- rownames(sce)
    colData(sce)$cell_type <- sce$cell_type
    saveRDS(sce, file="{sce}")

    png("{umap}")
    DimPlot(object = seurat, reduction = "umap", group.by = "cell_type")
    dev.off()

    png("{tsne}")
    DimPlot(object = seurat, reduction = "tsne", group.by = "cell_type")
    dev.off()

    """
    integrate_script = os.path.join(".cache/integration_seurat.R")
    output = open(integrate_script,"w")
    output.write(rcode.format(object_list=",".join(object_list), rdata=rdata, sce=sce_cached,umap=umap,tsne=tsne))
    output.close()
    subprocess.call(["Rscript","{}".format(integrate_script)])
    shutil.copyfile(rdata, integrated_seurat)
    shutil.copyfile(sce_cached, integrated_sce)
    shutil.copyfile(tsne, integrated_tsne)
    shutil.copyfile(umap, integrated_umap)

def RunHarmonyIntegration(sces, integrated_harmony, integrated_sce, integrated_tsne, integrated_umap, merged):
    rdata = os.path.join(config.jobpath,"results","integrated_harmony_seurat.rdata")
    sce_cached = os.path.join(config.jobpath,"results","integrated_harmony_sce.rdata")
    umap = os.path.join(config.jobpath,"results","integrated_harmony_tsne.rdata")
    tsne = os.path.join(config.jobpath,"results","integrated_harmony_umap.rdata")
    merged = os.path.join(config.jobpath,"results","sce_merged.rdata")
    object_list = []
    rcode = """
    library(Seurat)
    library(SingleCellExperiment)
    library(harmony)
    library(scater)
    """
    for idx, object in sces.items():
        seurat_obj = "seurat{}".format(idx)
        object_list.append(seurat_obj)
        load = """
        sce <- readRDS("{object}")
        colnames(sce) <- colData(sce)$Barcode
        sce$sample <- "{idx}"
        rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ensembl_gene_id, rownames(sce))
        seurat{idx} <- as.Seurat(sce, counts = "counts", data = "logcounts")
        """.format(idx=idx,object=object)
        rcode += load
    rcode += """
    merged <- merge(seurat1, y = c({object_list}), project = "pipeline_run")
    saveRDS(merged, file="{merged}")
    merged <- NormalizeData(merged)
    merged <- FindVariableFeatures(merged)
    merged <- ScaleData(merged)
    merged <- RunPCA(merged, verbose = FALSE)
    integrated <- RunHarmony(merged, c("sample","cell_type"), theta=c(0.1,1), lambda=c(4,1))
    integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:50)
    saveRDS(integrated, file="{rdata}")
    sce <- as.SingleCellExperiment(integrated)
    saveRDS(sce, file="{sce_cached}.rdata")
    """
    integrate_script = os.path.join(".cache/integration_harmony.R")
    output = open(integrate_script,"w")
    output.write(rcode.format(object_list=",".join(object_list), rdata=rdata, sce_cached=sce_cached,umap=umap,tsne=tsne,merged=merged))
    output.close()
    cmd = """Rscript {script}""".format(script=integrate_script)
    subprocess.call(cmd.split())
    shutil.copyfile(rdata, integrated_harmony)
    shutil.copyfile(sce_cached, integrated_sce)
    shutil.copyfile(tsne, integrated_tsne)
    shutil.copyfile(umap, integrated_umap)
    shutil.copyfile(merged, sce_merged)

def RunScanoramaIntegration(merged, integrated_sce, integrated_tsne, integrated_umap):
    rdata = os.path.join(config.jobpath,"results","integrated_scanorama_sce.rdata")
    umap = os.path.join(config.jobpath,"results","integrated_scanorama_tsne.rdata")
    tsne = os.path.join(config.jobpath,"results","integrated_scanorama_umap.rdata")
    batch_correct_method = """
    batch_correct <- function(sce, batch_col, method = "scanorama") {
      ## Feature selection
      batches <- unique(colData(sce)[,batch_col])

      fit <- trendVar(sce, parametric=TRUE, use.spikes = FALSE)
      decomp <- decomposeVar(sce, fit)
      decomp$Symbol <- rowData(sce)$Symbol

      decomp_chosen <- decomp %>% subset(bio > 0)
      chosen <- rownames(decomp_chosen)

      chosen_features <- rep(list(chosen), length(batches))

      ## Batch correction
      if (method == "scanorama") {
        norm_matrices <- lapply(batches, function(bat) {
          sce_sub <- sce[,colData(sce)[,batch_col] == bat]
          norm_data <- t(as.matrix(logcounts(sce_sub)))
          return(norm_data)
        })

        indexes <- lapply(batches, function(bat) {
          which(colData(sce)[,batch_col] == bat)
        })

        scanorama <- reticulate::import('scanorama')
        print("imported")
        # Integration and batch correction
        integrated_corrected_data <- scanorama$correct(norm_matrices,
                                                       chosen_features,
                                                       return_dimred = TRUE,
                                                       return_dense = TRUE)
        print("Finished")
        scanorama_mat <- matrix(NA, nrow = sum(sapply(indexes, length)), ncol = 100)
        scanorama_mat_expr <- matrix(NA, nrow = sum(sapply(indexes, length)), ncol = ncol(integrated_corrected_data[[2]][[1]]))
        for (i in seq_along(indexes)) {
          idx <- indexes[[i]]
          scanorama_mat[idx,] <- integrated_corrected_data[[1]][[i]]
          scanorama_mat_expr[idx,] <- integrated_corrected_data[[2]][[i]]
        }
        reducedDim(sce, "scanorama_int") <- scanorama_mat
        print("reducedDim")
        sce_sel <- sce
        sce_sel <- sce_sel[integrated_corrected_data[[3]],]
        assay(sce_sel, "scanorama") <- t(scanorama_mat_expr)
        pca_res <- pca2(sce_sel, ntop = 1000, ncomponents = 50, exprs_values = "scanorama")
        reducedDim(sce_sel, "PCA2") <- pca_res$x

        sce_sel <- runTSNE(sce_sel, use_dimred = "PCA2", n_dimred = 50)
        sce_sel <- runUMAP(sce_sel, use_dimred = "PCA2", n_dimred = 50)
        print("select")
        reducedDim(sce, "scanorama_PCA") <- reducedDim(sce_sel, "PCA2")
        print("here222")
        reducedDim(sce, "scanorama_TSNE") <- reducedDim(sce_sel, "TSNE")
        reducedDim(sce, "scanorama_UMAP") <- reducedDim(sce_sel, "UMAP")
        print("here")
        sce@metadata$batchcor_pca_res <- pca_res
        sce@metadata$batchcor_genes <- chosen
        saveRDS(sce_sel,file="SPECTRUM-OV-009-OMENT-ALL.rdata")
        print("done")
      } else {
        stop("Other methods not implemented.")
      }

      return(sce)
    }
    """
    integrate_script = os.path.join(".cache/integration_scanorama.R")
    output = open(integrate_script,"w")
    rcode.write("""
    library(SingleCellExperiment)
    library(scater)
    library(data.table)
    library(ggrepel)
    library(grid)
    library(scran)
    library(yaml)
    library(scales)
    library(feather)
    library(limma)
    library(scrna.utils)
    library(scrna.sceutils)
    library(cellassign)
    library(cellassign.utils)
    library(stringr)
    library(reticulate)
    use_python("/usr/lib/python")
    """)
    script.write(batch_correct_method)
    script.write("""
    merged <- readRDS("{merged}")
    """.format(patient=idx))
    script.write("""
    sce <- tryCatch({batch_correct(merged,"sample")},error=function(e) {batch_correct(merged,"sample")})
    print("passed")
    """)
    script.write("""
    rowData(sce) <- rowData(merged)
    rownames(sce) <- rownames(merged)
    saveRDS(sce, file="{sce_cached}")""".format(patient=idx))
    script.close()
    cmd = """/admin/lsf/10.1/linux3.10-glibc2.17-x86_64/bin/bsub -K -J "scanorama" -R "rusage[mem=4]" -R "select[type==CentOS7]" -W 03:00 -n 16 -o output -e error singularity exec /work/ceglian/images/scrna-r-base.img Rscript {script}""".format(script=integrate_script)
    subprocess.call(cmd.split())
    shutil.copyfile(rdata, integrated_sce)
    shutil.copyfile(tsne, integrated_tsne)
    shutil.copyfile(umap, integrated_umap)


def RunCollection(workflow):

    all_samples = open(config.samples, "r").read().splitlines()

    workflow.transform (
        name = "parse_collection",
        func = RunParse,
        args = (
            all_samples,
            pypeliner.managed.TempOutputFile("sample_path.json","sample")
        )
    )
    workflow.transform (
        name = "run_qc",
        func = RunQC,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("sample_path.json","sample"),
            pypeliner.managed.TempOutputFile("sce_raw.rdata","sample"),
            pypeliner.managed.TempOutputFile("filtered.rdata","sample"),
        )
    )

    workflow.transform (
        name = "run_cellassign",
        func = RunCellAssign,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("sample_path.json","sample"),
            pypeliner.managed.TempInputFile("filtered.rdata","sample"),
            pypeliner.managed.TempOutputFile("sce.rdata","sample")
        )
    )


    workflow.transform (
        name = "run_convert",
        func = RunConvert,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("sample_path.json","sample"),
            pypeliner.managed.TempInputFile("sce.rdata","sample"),
            pypeliner.managed.TempOutputFile("seurat.rdata","sample"),
        )
    )

    workflow.transform (
        name = "run_seurat_workflow",
        func = RunSeuratWorkflow,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("sample_path.json","sample"),
            pypeliner.managed.TempInputFile("seurat.rdata","sample"),
            pypeliner.managed.TempOutputFile("seurat_qcd.rdata","sample"),
            pypeliner.managed.TempOutputFile("sce_qcd.rdata","sample"),
        )
    )

    workflow.transform (
        name = "visualize_sample",
        func = RunSeuratViz,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("sample_path.json","sample"),
            pypeliner.managed.TempInputFile("seurat_qcd.rdata","sample"),
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
            pypeliner.managed.TempInputFile("sample_path.json","sample"),
            pypeliner.managed.TempInputFile("seurat_qcd.rdata","sample"),
            pypeliner.managed.TempOutputFile("markers.csv","sample"),
        )
    )

    # workflow.transform (
    #     name = "seurat_integrate",
    #     func = RunSeuratIntegration,
    #     args = (
    #         pypeliner.managed.TempInputFile("seurat_qcd.rdata","sample"),
    #         pypeliner.managed.TempOutputFile("integrated_seurat_seurat.rdata"),
    #         pypeliner.managed.TempOutputFile("integrated_seurat_sce.rdata"),
    #         pypeliner.managed.TempOutputFile("integrated_seurat_tsne.png"),
    #         pypeliner.managed.TempOutputFile("integrated_seurat_umap.png"),
    #     )
    # )

    workflow.transform (
        name = "harmony_integrate",
        func = RunHarmonyIntegration,
        args = (
            pypeliner.managed.TempInputFile("seurat_qcd.rdata","sample"),
            pypeliner.managed.TempOutputFile("integrated_harmony_seurat.rdata"),
            pypeliner.managed.TempOutputFile("integrated_harmony_sce.rdata"),
            pypeliner.managed.TempOutputFile("integrated_harmony_tsne.png"),
            pypeliner.managed.TempOutputFile("integrated_harmony_umap.png"),
            pypeliner.managed.TempOutputFile("merged_sce.rdata"),
        )
    )

    workflow.transform (
        name = "scanorama_integrate",
        func = RunScanoramaIntegration,
        args = (
            pypeliner.managed.TempInputFile("merged_sce.rdata"),
            pypeliner.managed.TempOutputFile("integrated_scanorama_sce.rdata"),
            pypeliner.managed.TempOutputFile("integrated_scanorama_tsne.png"),
            pypeliner.managed.TempOutputFile("integrated_scanorama_umap.png"),
        )
    )

    return workflow
