import numpy
import pandas as pd
from interface.singlecellexperiment import SingleCellExperiment
import anndata

class Scanpy(object):

    def __init__(self, counts, logcounts):
        self.counts = counts
        self.logcounts = logcounts

    @classmethod
    def fromRData(scanpy_class, rdata):
       
        if not rdata: raise Exception

        sce = SingleCellExperiment.fromRData(rdata)
        
        embedding_names =["UMAP"]
        embeddings = {}

        for embedding in embedding_names:
            embeddings[embedding] = sce.getReducedDims(embedding)

        counts = sce.assays["counts"]
        log_counts = sce.assays["logcounts"]

        row_data = pd.DataFrame(sce.rowData)
        col_data = pd.DataFrame(sce.colData) 

        row_data = row_data.set_index("Symbol")
        col_data = col_data.set_index("Barcode")

        counts_annobj = anndata.AnnData(X = counts.T, obs = col_data, var = row_data )
        logcounts_annobj = anndata.AnnData(X = log_counts.T, obs = col_data, var = row_data, uns = embeddings )
        
        scanpy = Scanpy(counts_annobj, logcounts_annobj)
        return scanpy
    
 
        
    def get_counts(self):
        return self.counts.X
    
    def get_norm_counts(self):
        return self.norm_counts.X
    
    def get_genes(self):
        #counts and logcounts have same gene lists
        return self.counts.var["Symbol"].tolist()
    
    def to_sce(self):
        return None

    def get_cell_assignments(self):
        barcodes = self.counts.obs["Barcode"]
        cell_types = self.counts.obs["cell_type"]
        assignments = {}

        for cell_type, barcode in zip(cell_types, barcodes):
            if cell_type not in assignments:
                assignments[cell_type]=[barcode]
            else:
                assignments[cell_type].append(barcode)
        
        return assignments

    def get_UMAP(self):
        return self.logcounts.uns["UMAP"]

def test():
    scanpy1 = Scanpy.fromRData("/home/abramsd/hack/hackathon_scrna/data/hgsoc_cd45p.rdata")
    print(scanpy1.counts.obs.index)
    print(scanpy1.counts.var.index)

#test()

#scanpy1 = Scanpy.fromRData("/home/abramsd/hack/hackathon_scrna/data/hgsoc_cd45p.rdata")
#print(scanpy1.get_cell_assignments())
