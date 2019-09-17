import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os
import json
import shutil
import subprocess
from collections import defaultdict
from scipy import stats
import math
import numpy

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

def RunExtract(sample_to_path, rdata_path):
    sample = json.loads(open(sample_to_path,"r").read())
    sampleid, path = list(sample.items()).pop()
    tenx_analysis = TenxAnalysis(path)
    tenx_analysis.load()
    tenx_analysis.extract()
    qc = QualityControl(tenx_analysis, sampleid)
    if not os.path.exists(qc.sce):
        qc.run(mito=config.mito)
    shutil.copyfile(qc.sce, rdata_path)

def RunCellAssign(sce, annot_sce):
    filtered_sce = os.path.join(os.path.split(sce)[0],"sce_cas.rdata")
    _rho_csv = os.path.join(os.path.split(sce)[0],"rho_csv_sub.csv")
    _fit = os.path.join(os.path.split(sce)[0],"fit.rdata")
    sampleid = sce.split("/")[-2]
    if not os.path.exists(filtered_sce):
        CellAssign.run(sce, config.rho_matrix, _fit, rho_csv=_rho_csv, lsf=False)
    shutil.copyfile(filtered_sce, annot_sce)
    path = os.getcwd()
    shutil.copyfile(_fit, os.path.join(path,"fit.rdata"))


def RunModeCopyNumber(copy_number_data):
    output = open(copy_number_data,"w")
    cell_to_clone = open(config.copy_cell_clones,"r").read().splitlines()
    cell_to_clone.pop(0)
    cell_to_clone = dict([x.split("\t") for x in cell_to_clone])
    copy_number_data = open(config.filtered_cell_cn,"r").read().splitlines()
    header = copy_number_data.pop(0).split("\t")
    output.write("chr,start,end,copy_number,clone\n")
    copy_number_mapping = defaultdict(lambda : defaultdict(list))
    for row in copy_number_data:
            row = row.split("\t")
            chrm, start, stop, width = row[:4]
            stop = str(int(start) + int(width))
            cells = dict(zip(header[4:],row[4:]))
            #ploidy = stats.mode([math.ceil(int(cn)) for cn in cells.values()])
            #ploidy = ploidy.mode[0]
            #if ploidy == 0:
            #       ploidy = 2
            for cell, copy_number in cells.items():
                #normalized_cn = math.ceil(int(copy_number) / (int(ploidy) / 2))
                copy_number_mapping[chrm+":"+start+":"+stop][cell_to_clone[cell]].append(int(copy_number))

    for loc in copy_number_mapping:
            chrm, start, end = loc.split(":")
            for clone, copy_numbers in copy_number_mapping[loc].items():
                    if clone == "None": continue
                    mode_copy_number = stats.mode(copy_numbers).mode[0]
                    line = ["chr"+chrm, start, end, str(int(mode_copy_number)), clone]
                    print(",".join(line))
                    output.write(",".join(line)+"\n")
    output.close()

def RunCloneAlignInput(sce, copy_number_data, clone_sce, cnv_mat, raw_cnv):
    rcode = """
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(dplyr)
    library(tidyverse)
    library(org.Hs.eg.db)
    library(SingleCellExperiment)

    sce <- readRDS('{sce}')
    rownames(sce) <- rowData(sce)$ensembl_gene_id
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    g <- genes(txdb, single.strand.genes.only=FALSE)
    entrezgene_ensembl_map <- as.list(org.Hs.egENSEMBL)
    entrezgene_ensembl_map <- lapply(entrezgene_ensembl_map, `[`, 1)
    cnv_data <- read_csv(file='{cnv}')
    cnv_gr <- makeGRangesFromDataFrame(cnv_data, keep.extra.columns = TRUE)
    olaps <- findOverlaps(g, cnv_gr)
    df_gene <- data_frame(entrezgene = names(g)[queryHits(olaps)],
               copy_number = mcols(cnv_gr)$copy_number[subjectHits(olaps)],
               clone = mcols(cnv_gr)$clone[subjectHits(olaps)])
    df_gene <- dplyr::filter(df_gene, entrezgene %in% names(entrezgene_ensembl_map)) %>%
      dplyr::mutate(ensembl_gene_id = unlist(entrezgene_ensembl_map[entrezgene])) %>%
      dplyr::select(ensembl_gene_id, entrezgene, copy_number, clone) %>%
      drop_na()
    df_gene <- df_gene %>% group_by(ensembl_gene_id) %>% tally() %>%
      filter(n == length(unique(df_gene$clone))) %>%
      inner_join(df_gene) %>%
      dplyr::select(-n)
    df_gene_expanded <- spread(df_gene, clone, copy_number)
    cnv_mat <- dplyr::select(df_gene_expanded, -ensembl_gene_id, -entrezgene) %>% as.matrix()
    rownames(cnv_mat) <- df_gene_expanded$ensembl_gene_id
    write.csv(cnv_mat,file='{raw_cnv}')
    keep_gene <- rowMins(cnv_mat) <= 6 & rowVars(cnv_mat) > 0
    cnv_mat <- cnv_mat[keep_gene,]
    common_genes <- intersect(rownames(cnv_mat), rownames(sce))
    sce <- sce[common_genes,]
    cnv_mat <- cnv_mat[common_genes,]
    saveRDS(file='{clone_sce}', sce)
    saveRDS(file='{cnv_mat}', cnv_mat)
    """
    path = os.path.split(sce)[0]
    convert_script = os.path.join(path,"build_input.R")
    output = open(convert_script,"w")
    output.write(rcode.format(sce=sce,cnv=copy_number_data,clone_sce=clone_sce,cnv_mat=cnv_mat,raw_cnv=raw_cnv))
    output.close()
    subprocess.call(["Rscript","{}".format(convert_script)])
    shutil.copyfile(clone_sce, os.path.join(path,"clone_input.rdata"))
    shutil.copyfile(cnv_mat, os.path.join(path,"cnv_mat.rdata"))

def RunCloneAlign(clone_sce, cnv_mat, annotated_sce, cal_fit):
    cwd = os.getcwd()
    annotated_sce_cached = os.path.join(os.path.split(clone_sce)[0],"clone_annotated_cached.rdata")
    cal_fit_cached = os.path.join(os.path.split(clone_sce)[0],"cal_cached.rdata")
    qplot_cached = os.path.join(os.path.split(clone_sce)[0],"qplot_cached.png")
    rcode = """
    library(clonealign)
    library(ggplot2)
    library(SingleCellExperiment)
    sce <- readRDS('{clone_sce}')
    sce <- sce[,colData(sce)$cell_type=="Ovarian.cancer.cell"]
    cnv <- readRDS('{cnv_mat}')
    sce <- sce[rowSums(counts(sce))>0,]
    cal <- clonealign(sce, cnv,remove_outlying_genes=FALSE)
    sce$clone <- cal$clone
    saveRDS(sce,file='{annotated_sce}')
    saveRDS(cal, file='{cal_fit}')
    """
    path = os.path.split(clone_sce)[0]
    run_script = os.path.join(path,"run_clonealign.R")
    output = open(run_script,"w")
    output.write(rcode.format(clone_sce=clone_sce,cnv_mat=cnv_mat,annotated_sce=annotated_sce_cached,cal_fit=cal_fit_cached))
    output.close()
    cmd = "/admin/lsf/10.1/linux3.10-glibc2.17-x86_64/bin/bsub -K -J cellassign -R rusage[mem=1] -n 20 -We 40 -o out -e err /opt/common/CentOS_7/singularity/singularity-3.0.1/bin/singularity exec --bind /admin --bind /opt --bind {} /home/ceglian/images/scrna-pipeline-clonealign.img".format(cwd)
    if not os.path.exists(annotated_sce_cached) or not os.path.exists(cal_fit_cached):
        subprocess.call(["Rscript","{}".format(run_script)])
    shutil.copyfile(annotated_sce_cached, annotated_sce)
    shutil.copyfile(cal_fit_cached, cal_fit)
    path = os.getcwd()
    shutil.copyfile(cal_fit_cached, os.path.join(path,"cal.rdata"))

def RunSeuratViz(seurat, umap, umap_celltype, umap_clone):
    marker_list = GeneMarkerMatrix.read_yaml(config.rho_matrix)
    markers = ["'" + marker + "'" for marker in marker_list.genes]
    umap_plot = os.path.join(os.path.split(seurat)[0],"umap.png")
    umap_celltype_plot = os.path.join(os.path.split(seurat)[0],"umap_celltype.png")
    umap_clone_plot = os.path.join(os.path.split(seurat)[0],"umap_clone.png")
    rcode = """
    library(Seurat)
    library(ggplot2)
    seurat <- readRDS("{seurat}")

    png("{umap}")
    DimPlot(object = seurat, reduction = "umap")
    dev.off()

    png("{umap_celltype}")
    DimPlot(object = seurat, reduction = "umap", group.by = "cell_type")
    dev.off()

    png("{umap_clone}")
    DimPlot(object = seurat, reduction = "umap", group.by = "clone")
    dev.off()
    """
    path = os.path.split(seurat)[0]
    qc_script = os.path.join(path,"viz.R")
    output = open(qc_script,"w")
    output.write(rcode.format(seurat=seurat, umap=umap_plot, umap_celltype=umap_celltype_plot,umap_clone=umap_clone_plot))
    output.close()
    if not os.path.exists(umap_clone_plot):
        subprocess.call(["Rscript","{}".format(qc_script)])
    shutil.copyfile(umap_plot, umap)
    shutil.copyfile(umap_celltype_plot, umap_celltype)
    shutil.copyfile(umap_clone_plot, umap_clone)


def RunConvert(sce_cell, sce_clone, seurat):
    seurat_cached = os.path.join(os.path.split(sce_cell)[0],"seurat_raw.rdata")
    rcode = """
    library(Seurat)
    library(SingleCellExperiment)
    library(scater)
    clone_sce <- readRDS('{sce_cell}')
    cell_sce <- readRDS('{sce_clone}')


    colnames(cell_sce) <- colData(cell_sce)$Barcode
    colnames(clone_sce) <- colData(clone_sce)$Barcode

    barcodes <- intersect(colnames(cell_sce),colnames(clone_sce))
    clone_sce <- clone_sce[,barcodes]
    cell_sce <- cell_sce[,barcodes]

    clone_sce$cell_type <- cell_sce$cell_type
    sce <- clone_sce
    rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ensembl_gene_id, rownames(sce))
    seurat <- as.Seurat(sce, counts = "counts", data = "logcounts")
    saveRDS(seurat,file='{seurat}')
    """
    path = os.path.split(sce_cell)[0]
    convert_script = os.path.join(path,"convert.R")
    output = open(convert_script,"w")
    output.write(rcode.format(sce_clone=sce_clone, sce_cell=sce_cell ,seurat=seurat_cached))
    output.close()
    if not os.path.exists(seurat_cached):
        subprocess.call(["Rscript","{}".format(convert_script)])
    shutil.copyfile(seurat_cached, seurat)

def RunSeuratWorkflow(seurat, qcd_seurat, qcd_sce):
    seurat_cached = os.path.join(os.path.split(seurat)[0],"seuret_annot.rdata")
    sce_cached = os.path.join(os.path.split(seurat)[0],"sce_annot.rdata")
    rcode = """
    library(Seurat)
    library(sctransform)
    library(SingleCellExperiment)
    seurat <- readRDS("{seurat}")
    seurat <- SCTransform(object = seurat)
    seurat <- RunPCA(object = seurat)
    seurat <- FindNeighbors(object = seurat)
    seurat <- FindClusters(object = seurat)
    seurat <- RunUMAP(object = seurat, reduction = "pca", dims = 1:20)
    saveRDS(seurat, file = '{qcd_seurat}')
    sce <- as.SingleCellExperiment(seurat)
    rowData(sce)$Symbol <- rownames(sce)
    saveRDS(sce, file="{qcd_sce}")
    """
    path = os.path.split(seurat)[0]
    qc_script = os.path.join(path,"qc.R")
    output = open(qc_script,"w")
    output.write(rcode.format(seurat=seurat, qcd_seurat=seurat_cached, qcd_sce=sce_cached))
    output.close()
    if not os.path.exists(seurat_cached) or not os.path.exists(sce_cached):
        subprocess.call(["Rscript", "{}".format(qc_script)])
    shutil.copyfile(seurat_cached, qcd_seurat)
    shutil.copyfile(sce_cached, qcd_sce)

def RunIntegration(seurats, integrated_seurat, integrated_sce, flowsort="full"):
    rdata = os.path.join(os.path.split(integrated_seurat)[0],"integrate_seurat_cached_{}.rdata".format(flowsort))
    sce_cached = os.path.join(os.path.split(integrated_seurat)[0],"integrate_sce_cached_{}.rdata".format(flowsort))
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
    """
    integrate_script = os.path.join(".cache/integration_{}.R".format(flowsort))
    output = open(integrate_script,"w")
    output.write(rcode.format(object_list=",".join(object_list), rdata=rdata, sce=sce_cached))
    output.close()
    if not os.path.exists(sce_cached):
        subprocess.call(["Rscript","{}".format(integrate_script)])
    shutil.copyfile(rdata, integrated_seurat)
    shutil.copyfile(sce_cached, integrated_sce)


def RunFigures(sce, umap_cell,umap_clone, umap_sample):
    path = os.path.split(sce)[0]
    umap_cell_cached = os.path.join(os.path.split(sce)[0],"umap_celltype_cached.png")
    umap_clone_cached = os.path.join(os.path.split(sce)[0],"umap_clone_cached.png")
    umap_sample_cached = os.path.join(os.path.split(sce)[0],"umap_sample_cached.png")
    rcode = """
    library(SingleCellExperiment)
    library(scater)
    library(stringr)
    sce <- readRDS('{sce}')

    sce$Sample <- lapply(sce$Sample, function (x) str_replace(x,".cache/",""))
    sce$Sample <- lapply(sce$Sample, function (x) str_replace(x,"/filtered_feature_bc_matrix",""))
    sce$Sample <- as.character(sce$Sample)
    png('{umap_clone}',width=1000,height=1000)
    plotReducedDim(sce, use_dimred = "UMAP", colour_by = "clone") +
        xlab('UMAP-1') +
        ylab("UMAP-2") +
        guides(fill = guide_legend(title = "Cell Type")) +
        theme_bw() +     scale_fill_manual(values=c("#58d1eb","#000000","#df4173","#f4005f","gray","#58d1eb","#98e024","#000000","navy","#58d1eb","#000000","#df4173","#f4005f","gray","#58d1eb","#98e024","#000000","navy"))  + ggtitle("Clone") + theme(plot.title=element_text(size=19, face = "bold"),axis.title.x=element_text(size=15, face = "bold"),axis.title.y=element_text(size=15, face = "bold"), legend.text=element_text(size=15, face = "bold"),axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15),legend.title=element_text(size=15))
    dev.off()
    png('{umap_sample}',width=1000,height=1000)
    plotReducedDim(sce, use_dimred = "UMAP", colour_by = "Sample") +
        xlab('UMAP-1') +
        ylab("UMAP-2") +
        guides(fill = guide_legend(title = "Cell Type")) +
        theme_bw() + scale_fill_manual(values=c("#58d1eb","#000000","#df4173","#f4005f","gray","#58d1eb","#98e024","#000000","navy","#58d1eb","#000000","#df4173","#f4005f","gray","#58d1eb","#98e024","#000000","navy"))  + ggtitle("Sample") + theme(plot.title=element_text(size=19, face = "bold"),axis.title.x=element_text(size=15, face = "bold"),axis.title.y=element_text(size=15, face = "bold"), legend.text=element_text(size=15, face = "bold"),axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15),legend.title=element_text(size=15))
    dev.off()
    png('{umap_celltype}',width=1000,height=1000)
    plotReducedDim(sce, use_dimred = "UMAP", colour_by = "cell_type") +
        xlab('UMAP-1') +
        ylab("UMAP-2") +
        guides(fill = guide_legend(title = "Cell Type")) +
        theme_bw() + scale_fill_manual(values=c("#58d1eb","#000000","#df4173","#f4005f","gray","#58d1eb","#98e024","#000000","navy","#58d1eb","#000000","#df4173","#f4005f","gray","#58d1eb","#98e024","#000000","navy"))  + ggtitle("Cell Type") + theme(plot.title=element_text(size=19, face = "bold"),axis.title.x=element_text(size=15, face = "bold"),axis.title.y=element_text(size=15, face = "bold"), legend.text=element_text(size=15, face = "bold"),axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15),legend.title=element_text(size=15))
    dev.off()
    """
    run_script = os.path.join(path,"run_figures.R")
    output = open(run_script,"w")
    output.write(rcode.format(sce=sce,umap_clone=umap_clone_cached,umap_celltype=umap_cell_cached,umap_sample=umap_sample_cached))
    output.close()
    subprocess.call(["Rscript","{}".format(run_script)])
    shutil.copyfile(umap_cell_cached, umap_cell)
    shutil.copyfile(umap_clone_cached, umap_clone)
    shutil.copyfile(umap_sample_cached, umap_sample)
    path = os.getcwd()
    shutil.copyfile(sce, os.path.join(path,"sce.rdata"))
    shutil.copyfile(umap_clone_cached, os.path.join(path,"umap_clone.png"))
    shutil.copyfile(umap_cell_cached, os.path.join(path,"umap_cell.png"))
    shutil.copyfile(umap_sample_cached, os.path.join(path,"umap_sample.png"))

def RunCloneAlignWorkflow(workflow):
    print("Creating workflow.")
    all_samples = open(config.samples, "r").read().splitlines()
    all_samples = [sample.strip() for sample in all_samples if sample.strip() != ""]
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
        )
    )
    workflow.transform (
        name = "run_cellassign",
        func = RunCellAssign,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("sample.rdata","sample"),
            pypeliner.managed.TempOutputFile("cell_annotated.rdata","sample"),
        )
    )
    workflow.transform (
        name = "run_modecn",
        func = RunModeCopyNumber,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempOutputFile("copy_number_data.csv","sample"),
        )
    )
    workflow.transform (
        name = "run_clonealigninputs",
        func = RunCloneAlignInput,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("cell_annotated.rdata","sample"),
            pypeliner.managed.TempInputFile("copy_number_data.csv","sample"),
            pypeliner.managed.TempOutputFile("clone.rdata","sample"),
            pypeliner.managed.TempOutputFile("cnv.rdata","sample"),
            pypeliner.managed.TempOutputFile("rawcnv.rdata","sample"),
        )
    )
    workflow.transform (
        name = "run_clonealign",
        func = RunCloneAlign,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("clone.rdata","sample"),
            pypeliner.managed.TempInputFile("cnv.rdata","sample"),
            pypeliner.managed.TempOutputFile("clone_annotated.rdata","sample"),
            pypeliner.managed.TempOutputFile("cal.rdata","sample"),
        )
    )


    if len(all_samples) > 1:
        workflow.transform (
            name = "run_convert",
            func = RunConvert,
            axes = ('sample',),
            args = (
                pypeliner.managed.TempInputFile("clone_annotated.rdata","sample"),
                pypeliner.managed.TempInputFile("cell_annotated.rdata","sample"),
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
                pypeliner.managed.TempOutputFile("sce_qcd.rdata","sample"),
            )
        )

        workflow.transform (
            name = "visualize_sample",
            func = RunSeuratViz,
            axes = ('sample',),
            args = (
                pypeliner.managed.TempInputFile("seurat_qcd.rdata","sample"),
                pypeliner.managed.TempOutputFile("seurat_umap.png","sample"),
                pypeliner.managed.TempOutputFile("seurat_umap_celltype.png","sample"),
                pypeliner.managed.TempOutputFile("seurat_umap_clone.png","sample"),
            )
        )

        workflow.transform (
            name = "integrate",
            func = RunIntegration,
            args = (
                pypeliner.managed.TempInputFile("seurat_qcd.rdata","sample"),
                pypeliner.managed.TempOutputFile("seurat_integrated.rdata"),
                pypeliner.managed.TempOutputFile("sce_integrated.rdata"),
            )
        )

        workflow.transform (
            name = "run_figures",
            func = RunFigures,
            args = (
                pypeliner.managed.TempInputFile("sce_integrated.rdata"),
                pypeliner.managed.TempOutputFile("umap_cell.png"),
                pypeliner.managed.TempOutputFile("umap_clone.png"),
                pypeliner.managed.TempOutputFile("umap_sample.png"),
            )
        )
    else:
        workflow.transform (
            name = "run_figures_single_sample",
            func = RunFigures,
            axes = ('sample',),
            args = (
                pypeliner.managed.TempInputFile("clone_annotated.rdata","sample"),
                pypeliner.managed.TempOutputFile("umap_cell.png","sample"),
                pypeliner.managed.TempOutputFile("umap_clone.png","sample"),
                pypeliner.managed.TempOutputFile("umap_sample.png","sample"),
            )
        )
    return workflow
