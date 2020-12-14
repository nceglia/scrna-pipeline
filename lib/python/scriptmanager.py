import os

class Script(object):
    def __init__(self, name):
        self.name = name
        self.args = dict()
        self.outs = dict()

    def set_path(self, path):
        self.path = path
        ext = os.path.splitext(path)[1]
        if ext == ".py":
            self.interpreter = "python3"
        elif ext == ".R":
            self.interpreter = "Rscript"
        else:
            self.interpreter = "sh"

    def set_args(self, position, value):
        self.args[position] = value

    def set_outs(self, position, value):
        self.outs[position] = value

    def set_interpreter(self, interpreter):
        self.interpreter = interpreter

class ScriptManager(object):

    def __init__(self):
        self.current = os.path.split(os.path.realpath(__file__))[0]

    def qc(self):
        script = Script("qc")
        script.set_args(1, "matrix")
        script.set_args(2, "mito")
        script.set_args(3, "ncounts")
        script.set_outs(4, "seurat")
        script.set_outs(5, "sce")
        script.set_outs(6, "raw_sce")
        script.set_outs(7, "matrix")
        script.set_args(8, "nfeatures")
        script.set_path(os.path.join(self.current,"../R/qc.R"))
        return script

    def buildmatrix(self):
        script = Script("buildmatrix")
        script.set_args(1, "inventory")
        script.set_outs(2, "matrix")
        script.set_outs(3, "features")
        script.set_outs(4, "barcodes")
        script.set_outs(5, "project_matrix")
        script.set_path(os.path.join(self.current,"../python/buildmatrix.py"))
        return script

    def detectdoublets(self):
        script = Script("detectdoublets")
        script.set_args(1, "sce")
        script.set_outs(2, "csv")
        script.set_path(os.path.join(self.current,"../python/detectdoublets.py"))
        return script

    def filterdoublets(self):
        script = Script("filterdoublets")
        script.set_args(1, "seurat")
        script.set_args(2, "csv")
        script.set_args(3, "sample")
        script.set_args(4, "score")
        script.set_outs(5, "qcd_seurat")
        script.set_path(os.path.join(self.current,"../R/filterdoublets.R"))
        return script

    def cellcyclescore(self):
        script = Script("cellcyclescore")
        script.set_args(1, "qcd_seurat")
        script.set_outs(2, "qcd_scored_seurat")
        script.set_outs(3, "qcd_scored_sce")
        script.set_path(os.path.join(self.current,"../R/cellcycle.R"))
        return script

    def mergesamples(self):
        script = Script("mergesamples")
        script.set_outs(1, "merged_tsv")
        script.set_outs(2, "merged_seurat")
        script.set_args(3, "metadata")
        script.set_path(os.path.join(self.current,"../R/mergesamples.R"))
        return script

    def mergebatches(self):
        script = Script("mergebatches")
        script.set_outs(1, "batch_tsv")
        script.set_outs(2, "batch_merged")
        script.set_outs(3, "celltype_csv")
        script.set_args(4, "project_matrix")
        script.set_path(os.path.join(self.current,"../python/mergebatches.py"))
        return script

    def assigncelltypes(self):
        script = Script("assigncelltypes")
        script.set_args(1, "merged_seurat")
        script.set_args(2, "marker_csv")
        script.set_outs(3, "probabilities")
        script.set_outs(4, "annotated_seurat")
        script.set_outs(5, "batch_report")
        script.set_outs(6, "batch_csv")
        script.set_outs(7, "annotated_sce")
        script.set_path(os.path.join(self.current,"../R/cellassign.R"))
        return script

    def batchcorrection(self):
        script = Script("batchcorrection")
        script.set_args(1, "batch_merged")
        script.set_outs(2, "integrated_seurat")
        script.set_outs(3, "project_figure")
        script.set_outs(4, "ct_markers")
        script.set_outs(5, "cell_annotations")
        script.set_args(6, "celltype_csv")
        script.set_path(os.path.join(self.current,"../R/batchcorrection.R"))
        return script

    def subsetcelltype(self):
        script = Script("subsetcelltype")
        script.set_args(1, "integrated_seurat")
        script.set_args(2, "celltype")
        script.set_outs(3, "celltype_seurat")
        script.set_outs(4, "markers")
        script.set_outs(5, "markers_tsv")
        script.set_outs(6, "cells_tsv")
        script.set_args(7, "resolution")
        script.set_path(os.path.join(self.current,"../R/subsetcelltype.R"))
        return script

    def celltypemarkers(self):
        script = Script("celltypemarkers")
        script.set_args(1, "yaml")
        script.set_outs(2, "marker_csv")
        script.set_path(os.path.join(self.current,"../python/parseyaml.py"))
        return script

    def summarizeqc(self):
        script = Script("qcreport")
        script.set_args(1, "seurat")
        script.set_args(2, "sce")
        script.set_args(3, "csv")
        script.set_args(4, "markers")
        script.set_args(5, "sample")
        script.set_args(6, "mito")
        script.set_args(7, "score")
        script.set_outs(8, "qc_report")
        script.set_path(os.path.join(self.current,"../R/summarizeqc.R"))
        return script
