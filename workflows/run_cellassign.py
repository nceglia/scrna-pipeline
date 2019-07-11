import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os
import pickle
import shutil

from interface.tenxanalysis import TenxAnalysis
from software.cellassign import CellAssign
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl
from interface.genemarkermatrix import GeneMarkerMatrix
from utils.plotting import celltypes, tsne_by_cell_type, umap_by_cell_type
from software.dropletutils import DropletUtils

from utils.config import Configuration

config = Configuration()

def Run(sampleid, raw_sce, sce_cas):
    if not os.path.exists(".cache/{}/celltypes.rdata".format(sampleid)):
        CellAssign.run(raw_sce, config.rho_matrix, ".cache/{}/celltypes.rdata".format(sampleid))
    shutil.copyfile(".cache/{}/sce_cas.rdata".format(sampleid), sce_cas)

def Analysis(sampleid, sce_cas, celltypes, tsne, umap):
    filtered_sce = sce_cas
    cellassign_analysis = ".cache/{}/cellassignanalysis/".format(sampleid)
    if not os.path.exists(cellassign_analysis):
        os.makedirs(cellassign_analysis)
    pyfit = os.path.join(".cache/{}/cell_types.pkl".format(sampleid))
    assert os.path.exists(pyfit), "No Pyfit Found."
    pyfit = pickle.load(open(pyfit,"rb"))
    marker_list = GeneMarkerMatrix.read_yaml(config.rho_matrix)
    cell_types = marker_list.celltypes()
    if "B cell" not in cell_types: cell_types.append("B cell")
    celltypes(pyfit, sampleid, cellassign_analysis, known_types=cell_types)
    tsne_by_cell_type(filtered_sce, pyfit, sampleid, cellassign_analysis, known_types=cell_types)
    umap_by_cell_type(filtered_sce, pyfit, sampleid, cellassign_analysis, known_types=cell_types)
    _celltypes = os.path.join(cellassign_analysis, "cell_types.png")
    _tsne = os.path.join(cellassign_analysis, "tsne_by_cell_type.png")
    _umap = os.path.join(cellassign_analysis, "umap_by_cell_type.png")
    shutil.copyfile(_celltypes, celltypes)
    shutil.copyfile(_umap, umap)
    shutil.copyfile(_tsne, tsne)

def RunCellAssign(sampleid, workflow):
    workflow.transform (
        name = "cellassign",
        func = Run,
        args = (
            sampleid,
            pypeliner.managed.TempInputFile("raw_sce.rdata"),
            pypeliner.managed.TempOutputFile("sce_cas.rdata")
        )
    )
    workflow.transform (
        name = "cellassignanalysis",
        func = Analysis,
        args = (
            sampleid,
            pypeliner.managed.TempInputFile("sce_cas.rdata"),
            pypeliner.managed.TempOutputFile("celltypes.png"),
            pypeliner.managed.TempOutputFile("tsne_by_celltype.png"),
            pypeliner.managed.TempOutputFile("umap_by_celltype.png"),
        )
    )
    return workflow

if __name__ == '__main__':
    sampleid = sys.argv[1]
    Run(sampleid, "qc.complete", "cellassign.complete")
