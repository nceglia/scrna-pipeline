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

from interface.singlecellexperiment import SingleCellExperiment

from utils.config import Configuration, write_config

config = Configuration()


def RunSeuratIntegration(sample_paths, integrated_seurat, integrated_sce, integrated_umap):
    if not os.path.exists(".cache"):
        os.makedirs(".cache")
    if not os.path.exists('results'):
        os.makedirs("results")
    rdata = os.path.join(config.jobpath,"results","integrated_seurat_seurat.rdata")
    sce_cached = os.path.join(config.jobpath,"results","integrated_seurat_sce.rdata")
    umap = os.path.join(config.jobpath,"results","integrated_seurat_umap.rdata")
    object_list = []
    rcode = """
    library(Seurat)
    library(SingleCellExperiment)
    """
    for idx, object in sample_paths.items():
        seurat_obj = "seurat{}".format(idx)
        object_list.append(seurat_obj)
        load = """
    sce{id} <- readRDS("{object}")
    seurat{id} <- as.Seurat(sce{id}, counts = "counts", data = "logcounts")
    seurat{id} <- SCTransform(seurat{id})
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

    """
    integrate_script = os.path.join(".cache/integration_seurat.R")
    output = open(integrate_script,"w")
    output.write(rcode.format(object_list=",".join(object_list), rdata=rdata, sce=sce_cached,umap=umap))
    output.close()
    subprocess.call(["Rscript","{}".format(integrate_script)])
    shutil.copyfile(rdata, integrated_seurat)
    shutil.copyfile(sce_cached, integrated_sce)
    shutil.copyfile(umap, integrated_umap)

def RunHarmonyIntegration(sample_paths, integrated_harmony, integrated_sce, integrated_umap, merged):
    if not os.path.exists(".cache"):
        os.makedirs(".cache")
    if not os.path.exists('results'):
        os.makedirs("results")
    rdata = os.path.join(config.jobpath,"results","integrated_harmony_seurat.rdata")
    sce_cached = os.path.join(config.jobpath,"results","integrated_harmony_sce.rdata")
    umap = os.path.join(config.jobpath,"results","integrated_harmony_umap.png")
    merged = os.path.join(config.jobpath,"results","sce_merged.rdata")
    object_list = []
    rcode = """
    library(Seurat)
    library(SingleCellExperiment)
    library(harmony)
    library(scater)
    """
    for idx, object in sample_paths.items():
        seurat_obj = "seurat{}".format(idx)
        object_list.append(seurat_obj)
        load = """
        sce <- readRDS("{object}")
        colnames(sce) <- colData(sce)$Barcode
        sce$sample <- "{idx}"
        seurat{idx} <- as.Seurat(sce, counts = "counts", data = "logcounts")
        """.format(idx=idx,object=object)
        rcode += load
    init_object = object_list.pop(0)
    rcode += """
    merged <- merge({init_object}, y = c({object_list}), project = "pipeline_run")
    saveRDS(merged, file="{merged}")
    merged <- NormalizeData(merged)
    merged <- FindVariableFeatures(merged)
    merged <- ScaleData(merged)
    merged <- RunPCA(merged, verbose = FALSE)
    integrated <- RunHarmony(merged, c("sample","cell_type"), theta=c(0.1,1), lambda=c(4,1))
    integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:50)
    saveRDS(integrated, file="{rdata}")
    sce <- as.SingleCellExperiment(integrated)
    saveRDS(sce, file="{sce_cached}")

    png("{umap}")
    DimPlot(object = seurat, reduction = "umap", group.by = "cell_type")
    dev.off()
    """
    integrate_script = os.path.join(".cache/integration_harmony.R")
    output = open(integrate_script,"w")
    output.write(rcode.format(init_object=init_object,object_list=",".join(object_list), rdata=rdata, sce_cached=sce_cached,umap=umap,merged=merged))
    output.close()
    cmd = """Rscript {script}""".format(script=integrate_script)
    subprocess.call(cmd.split())
    shutil.copyfile(rdata, integrated_harmony)
    shutil.copyfile(sce_cached, integrated_sce)
    shutil.copyfile(umap, integrated_umap)
    shutil.copyfile(merged, sce_merged)

def RunScanoramaIntegration(merged, integrated_sce, integrated_umap):
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
    print(config.samples)
    all_samples = json.loads(open(config.samples, "r").read())

    workflow.transform (
        name = "seurat_integrate",
        func = RunSeuratIntegration,
        args = (
            all_samples,
            pypeliner.managed.TempOutputFile("integrated_seurat_seurat.rdata"),
            pypeliner.managed.TempOutputFile("integrated_seurat_sce.rdata"),
            pypeliner.managed.TempOutputFile("integrated_seurat_umap.png"),
        )
    )

    workflow.transform (
        name = "harmony_integrate",
        func = RunHarmonyIntegration,
        args = (
            all_samples,
            pypeliner.managed.TempOutputFile("integrated_harmony_seurat.rdata"),
            pypeliner.managed.TempOutputFile("integrated_harmony_sce.rdata"),
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
            pypeliner.managed.TempOutputFile("integrated_scanorama_umap.png"),
        )
    )

    return workflow
