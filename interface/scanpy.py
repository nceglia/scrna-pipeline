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
       
        if not rdata: raise Exception

        sce = SingleCellExperiment.fromRData(rdata)
        print( sce.assays.keys())
        counts = sce.assays["counts"]
        log_counts = sce.assays["logcounts"]

        row_data = pd.DataFrame(sce.rowData)
        col_data = pd.DataFrame(sce.colData)

        counts_annobj = anndata.AnnData(X = counts, obs = row_data, var = col_data )
        logcounts_annobj = anndata.AnnData(X = log_counts, obs = row_data, var = col_data )
        scanpy = Scanpy(annobjannobj, logcounts_annobj)
        return scanpy

    def get_counts():
        return self.counts.X
    
    def get_norm_counts():
        return self.norm_counts.X
    

    def to_sce(self):
        return None

scanpy1 = Scanpy.fromRData("/home/abramsd/hack/hackathon_scrna/data/lung_1.rdata")
scanpy1.get_counts()

