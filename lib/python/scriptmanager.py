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
        script.set_outs(3, "seurat")
        script.set_outs(4, "sce")
        script.set_outs(5, "raw_sce")
        script.set_path(os.path.join(self.current,"../R/qc.R"))
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
        script.set_outs(3, "qcd_scored_seurat")
        script.set_path(os.path.join(self.current,"../R/cellcycle.R"))
        return script

    def mergesamples(self):
        script = Script("mergesamples")
        script.set_outs(1, "merged_tsv")
        script.set_outs(2, "merged_seurat")
        script.set_path(os.path.join(self.current,"../R/mergesamples.R"))
        return script

    def mergebatches(self):
        script = Script("mergebatches")
        script.set_outs(1, "batch_tsv")
        script.set_outs(2, "batch_merged_seurat")
        script.set_outs(3, "celltype_csv")
        script.set_path(os.path.join(self.current,"../R/mergebatches.R"))
        return script

    def assigncelltypes(self):
        script = Script("assigncelltypes")
        script.set_args(1, "merged_seurat")
        script.set_args(2, "marker_csv")
        script.set_outs(3, "probabilities")
        script.set_outs(4, "annotated_seurat")
        script.set_outs(5, "batch_report")
        script.set_outs(6, "batch_csv")
        script.set_path(os.path.join(self.current,"../R/cellassign.R"))
        return script

    def batchcorrection(self):
        script = Script("batchcorrection")
        script.set_args(1, "batch_merged_seurat")
        script.set_outs(2, "integrated_seurat")
        script.set_outs(3, "project_figure")
        script.set_outs(4, "ct_markers")
        script.set_path(os.path.join(self.current,"../R/batchcorrection.R"))
        return script

    def enrichmentnetwork(self):
        script = Script("enrichmentnetwork")
        script.set_args(1, "ct_markers")
        script.set_args(2, "celltype")
        script.set_outs(3, "pathway_network")
        script.set_args(4, "gmt")
        script.set_outs(5, "enriched_pathways")
        script.set_path(os.path.join(self.current,"../python/pathwayanalysis.py"))
        return script

    def subsetcelltype(self):
        script = Script("subsetcelltype")
        script.set_args(1, "batch_merged_seurat")
        script.set_args(2, "celltype")
        script.set_outs(3, "celltype_seurat")
        script.set_outs(4, "cell_umap")
        script.set_outs(5, "markers")
        script.set_path(os.path.join(self.current,"../R/subsetcelltype.R"))
        return script

    def cellcellinteractions(self):
        script = Script("cellcellinteractions")
        script.set_args(1, "sender")
        script.set_args(2, "receiver")
        script.set_args(3, "sender_obj")
        script.set_args(4, "receiver_obj")
        script.set_args(5, "geneset")
        script.set_args(6, "pathway")
        script.set_outs(7, "interactions")
        script.set_path(os.path.join(self.current,"../R/cellcellinteractions.R"))
        return script


    def subtypescores(self):
        script = Script("subtypescores")
        script.set_args(1, "celltype_obj")
        script.set_args(2, "subtype")
        script.set_args(3, "geneset")
        script.set_outs(7, "subtype_scores")
        script.set_path(os.path.join(self.current,"../R/subtypescores.R"))
        return script

    def celltypemarkers(self):
        script = Script("celltypemarkers")
        script.set_args(1, "yaml")
        script.set_outs(2, "marker_csv")
        script.set_path(os.path.join(self.current,"../python/parseyaml.py"))
        return script

    def generatereport(self):
        script = Script("generatereport")
        script.set_args(1, "qc_report")
        script.set_args(2, "project")
        script.set_args(3, "batch_report")
        script.set_args(4, "project_figure")
        script.set_args(8, "celltype_umaps")
        script.set_args(9, "interactions")
        script.set_args(10, "networks")
        script.set_args(11, "markers")
        script.set_args(12, "batch_markers")
        script.set_args(13, "subtype_scores")
        script.set_args(14, "enriched_pathways")
        script.set_args(15, "ct_markers")
        script.set_outs(16, "report")
        script.set_path(os.path.join(self.current,"../python/report.py"))
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
