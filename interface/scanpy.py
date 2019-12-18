import numpy
import pandas as pd
from interface.singlecellexperiment import SingleCellExperiment
import anndata

class Scanpy(object):

    def __init__(self, counts, norm_counts):
        self.counts = counts
        self.norm_counts = norm_counts

    @classmethod
    def fromRData(scanpy_class, rdata):
       
        sce = SingleCellExperiment.fromRData(rdata)
       
        counts = sce.assays["counts"]
        row_data = pd.DataFrame(sce.rowData)
        col_data = pd.DataFrame(sce.colData)
        annobj = anndata.AnnData(X = counts, obs = row_data, var = col_data )

        scanpy = Scanpy(annobj)
        return scanpy

    def get_counts  
    def to_sce(self):
        return None

scanpy1 = Scanpy.fromRData("/home/abramsd/hack/hackathon_scrna/data/lung_1.rdata")


