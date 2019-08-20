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
        CellAssign.run(sce, config.rho_matrix, _fit, rho_csv=_rho_csv)
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

def RunCloneAlignInput(sce, copy_number_data, clone_sce, cnv_mat):
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
    output.write(rcode.format(sce=sce,cnv=copy_number_data,clone_sce=clone_sce,cnv_mat=cnv_mat))
    output.close()
    subprocess.call(["Rscript","{}".format(convert_script)])
    shutil.copyfile(clone_sce, os.path.join(path,"clone_input.rdata"))
    shutil.copyfile(cnv_mat, os.path.join(path,"cnv_mat.rdata"))

def RunCloneAlign(clone_sce, cnv_mat, annotated_sce, cal_fit, qplot):
    annotated_sce_cached = os.path.join(os.path.split(clone_sce)[0],"clone_annotated_cached.rdata")
    cal_fit_cached = os.path.join(os.path.split(clone_sce)[0],"cal_cached.rdata")
    qplot_cached = os.path.join(os.path.split(clone_sce)[0],"qplot_cached.rdata")
    rcode = """
    library(clonealign)
    library(ggplot2)
    sce <- readRDS('{clone_sce}')
    cnv <- readRDS('{cnv_mat}')
    cal <- clonealign(sce, cnv,remove_outlying_genes=FALSE)
    sce$clone <- cal$clone
    saveRDS(sce,file='{annotated_sce}')
    saveRDS(cal, file='{cal_fit}')
    png('{qplot}')
    qplot(seq_along(cal$elbo), cal$elbo, geom = c("point", "line")) + labs(x = "Iteration", y = "ELBO")
    dev.off(0)
    """
    path = os.path.split(clone_sce)[0]
    run_script = os.path.join(path,"run_clonealign.R")
    output = open(run_script,"w")
    output.write(rcode.format(clone_sce=clone_sce,cnv_mat=cnv_mat,annotated_sce=annotated_sce_cached,cal_fit=cal_fit_cached, qplot=qplot_cached))
    output.close()
    if not os.path.exists(annotated_sce_cached) or not os.path.exists(cal_fit_cached):
        subprocess.call(["Rscript","{}".format(run_script)])
    shutil.copyfile(annotated_sce_cached, annotated_sce)
    shutil.copyfile(cal_fit_cached, cal_fit)
    shutil.copyfile(qplot_cached, qplot)
    path = os.getcwd()
    shutil.copyfile(cal_fit_cached, os.path.join(path,"cal.rdata"))
    shutil.copyfile(qplot_cached, os.path.join(path,"qplot.png"))

def RunEvaluation(annotated_sce, cal_fit, cnv_mat, evaluate_png):
    evaluate_png_cached = os.path.join(os.path.split(annotated_sce)[0],"evaluation_cached.png")
    rcode = """
    library(tidyverse)
    library(scater)
    library(data.table)
    library(broom)
    library(clonealign)
    sce <- readRDS('{annotated_sce}')
    cal <- readRDS('{cal_fit}')
    cnv <- readRDS('{cnv_mat}')""".format(annotated_sce=annotated_sce,cnv_mat=cnv_mat,cal_fit=cal_fit)

    rcode += """
    recompute_clone_assignment <- function(ca, clone_assignment_probability = 0.95) {
      clone_names <- colnames(ca$ml_params$clone_probs)
      clones <- apply(ca$ml_params$clone_probs, 1, function(r) {
        if(max(r) < clone_assignment_probability) {
          return("unassigned")
        }
        return(clone_names[which.max(r)])
      })
      ca$clone <- clones
      ca
    }

    ca <- recompute_clone_assignment(cal, 0.5)
    cnv <- dplyr::filter(cnv, !use_gene) %>%
      dplyr::rename(clone = cluster,
                    copy_number=median_cnmode) %>%
      dplyr::select(ensembl_gene_id, clone, copy_number) %>%
      spread(clone, copy_number)
    inferred_clones <- unique(ca$clone)
    inferred_clones <- setdiff(inferred_clones, "unassigned")

    collapsed_clones <- grepl("_", inferred_clones)

    if(any(collapsed_clones)) {
      for(i in which(collapsed_clones)) {
        cclone <- inferred_clones[i]
        uclones <- unlist(strsplit(cclone, "_"))
        new_clone <- rowMedians(cnv_mat[, uclones])
        cnv_mat <- cbind(cnv_mat, new_clone)
        colnames(cnv_mat)[ncol(cnv_mat)] <- cclone
      }
    }
    cnv_mat <- cnv_mat[, inferred_clones]
    cnv_mat <- cnv_mat[matrixStats::rowVars(cnv_mat) > 0,]
    cnv_mat <- cnv_mat[matrixStats::rowMaxs(cnv_mat) < 6,]

    sce <- sce[,ca$clone_fit$Barcode]

    sce <- sce[rowSums(as.matrix(counts(sce))) > 100, ]

    common_genes <- intersect(rownames(sce), rownames(cnv_mat))

    sce <- sce[common_genes,]
    cnv_mat <- cnv_mat[common_genes,]

    clones <- ca$clone

    assigned_cells <- clones != "unassigned"

    sce <- sce[, assigned_cells]
    clones <- clones[assigned_cells]

    logcs <- logcounts(sce)
    cnv_mat_full <- cnv_mat[, clones]

    test_estimates <- lapply(seq_len(nrow(sce)), function(i) {
      lc <- logcs[i,]
      cnv_dist <- cnv_mat_full[i,]
      tidy(lm(lc ~ cnv_dist))[2,]
    }) %>%
      bind_rows()


    cnv_mat_full <- cnv_mat[, sample(clones)]
    null_estimates <- lapply(seq_len(nrow(sce)), function(i) {
      lc <- logcs[i,]
      cnv_dist <- cnv_mat_full[i,]
      tidy(lm(lc ~ cnv_dist))[2,]
    }) %>%
      bind_rows()

    df <- bind_rows(
      dplyr::mutate(test_estimates, dist = "observed"),
      dplyr::mutate(null_estimates, dist = "null")
    )

    tt <- t.test(test_estimates$estimate, null_estimates$estimate)

    round2 <- function(x) format(round(x, 2), nsmall = 2)
    """
    rcode += """
    png('{evaluate_png}')
    ggplot(df, aes(x = dist, y = estimate)) +
      geom_boxplot(outlier.shape = NA, size = .4) +
      labs(x = "Distribution", y = "Coefficient expression ~ copy number",
           title = sample,
           subtitle = paste0("Genes not used by clonealign, p = ", round2(tt$p.value))) +
      theme_bw() +
      ylim(-.2, .2)
    dev.off()
    """.format(evaluate_png=evaluate_png_cached)
    path = os.path.split(annotated_sce)[0]
    run_script = os.path.join(path,"run_evaluation.R")
    output = open(run_script,"w")
    output.write(rcode)
    output.close()
    subprocess.call(["Rscript","{}".format(run_script)])
    shutil.copyfile(evaluate_png_cached, evaluate_png)


def RunFigures(clone_sce, cell_sce, result_sce, tsne, umap):
    path = os.path.split(cell_sce)[0]
    result_sce_cached = os.path.join(os.path.split(clone_sce)[0],"sce_cached.rdata")
    tsne_cached = os.path.join(os.path.split(clone_sce)[0],"tsne_cached.png")
    umap_cached = os.path.join(os.path.split(clone_sce)[0],"umap_cached.png")
    rcode = """
    library(SingleCellExperiment)
    library(scater)
    clone_sce <- readRDS('{clone_sce}')
    cell_sce <- readRDS('{cell_sce}')

    cell_sce <- cell_sce[,cell_sce$cell_type=="Ovarian.cancer.cell"]

    colnames(cell_sce) <- colData(cell_sce)$Barcode
    colnames(clone_sce) <- colData(clone_sce)$Barcode

    barcodes <- intersect(colnames(cell_sce),colnames(clone_sce))
    clone_sce <- clone_sce[,barcodes]
    cell_sce <- cell_sce[,barcodes]

    cell_sce$clone <- clone_sce$clone

    saveRDS(cell_sce,file='{result_sce}')

    png('{umap}')
    plotReducedDim(cell_sce, use_dimred = "UMAP", colour_by = "clone") +
        xlab('UMAP-1') +
        ylab("UMAP-2") +
        guides(fill = guide_legend(title = "Cell Type")) +
        theme_bw() + scale_fill_manual(values=c("#58d1eb","#000000","#df4173","#f4005f","gray","#58d1eb","#98e024","#000000","navy"))  + ggtitle("CloneAlign") + theme(plot.title=element_text(size=19, face = "bold"),axis.title.x=element_text(size=15, face = "bold"),axis.title.y=element_text(size=15, face = "bold"), legend.text=element_text(size=15, face = "bold"),axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15),legend.title=element_text(size=15))
    dev.off()

    png('{tsne}')
    plotReducedDim(cell_sce, use_dimred = "TSNE", colour_by = "clone") +
        xlab('UMAP-1') +
        ylab("UMAP-2") +
        guides(fill = guide_legend(title = "Cell Type")) +
        theme_bw() + scale_fill_manual(values=c("#58d1eb","#000000","#df4173","#f4005f","gray","#58d1eb","#98e024","#000000","navy"))  + ggtitle("CloneAlign") + theme(plot.title=element_text(size=19, face = "bold"),axis.title.x=element_text(size=15, face = "bold"),axis.title.y=element_text(size=15, face = "bold"), legend.text=element_text(size=15, face = "bold"),axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15),legend.title=element_text(size=15))
    dev.off()

    """
    run_script = os.path.join(path,"run_figures.R")
    output = open(run_script,"w")
    output.write(rcode.format(clone_sce=clone_sce,cell_sce=cell_sce,result_sce=result_sce_cached,tsne=tsne_cached,umap=umap_cached))
    output.close()
    subprocess.call(["Rscript","{}".format(run_script)])
    shutil.copyfile(result_sce_cached, result_sce)
    shutil.copyfile(umap_cached, umap)
    shutil.copyfile(tsne_cached, tsne)
    path = os.getcwd()
    shutil.copyfile(result_sce_cached, os.path.join(path,"sce.rdata"))
    shutil.copyfile(umap_cached, os.path.join(path,"umap.png"))
    shutil.copyfile(tsne_cached, os.path.join(path,"tsne.png"))

def RunCloneAlignWorkflow(workflow):
    print("Creating workflow.")
    all_samples = [config.prefix]
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
            pypeliner.managed.TempInputFile("sample.rdata","sample"),
            pypeliner.managed.TempInputFile("copy_number_data.csv","sample"),
            pypeliner.managed.TempOutputFile("clone.rdata","sample"),
            pypeliner.managed.TempOutputFile("cnv.rdata","sample"),
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
            pypeliner.managed.TempOutputFile("qplot.png","sample"),
        )
    )
    workflow.transform (
        name = "run_cloneeval",
        func = RunEvaluation,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("clone_annotated.rdata","sample"),
            pypeliner.managed.TempInputFile("cal.rdata","sample"),
            pypeliner.managed.TempInputFile("cnv.rdata","sample"),
            pypeliner.managed.TempOutputFile("clone_evaluation.png","sample"),
        )
    )
    workflow.transform (
        name = "run_figures",
        func = RunFigures,
        axes = ('sample',),
        args = (
            pypeliner.managed.TempInputFile("clone_annotated.rdata","sample"),
            pypeliner.managed.TempInputFile("cell_annotated.rdata","sample"),
            pypeliner.managed.TempOutputFile("sce.rdata","sample"),
            pypeliner.managed.TempOutputFile("tsne.png","sample"),
            pypeliner.managed.TempOutputFile("umap.png","sample"),
        )
    )
    return workflow
