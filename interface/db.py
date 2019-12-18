import os
import pandas as pd
import numpy


class Panglao(object):

    def __init__(self, file):
        csv = pd.read_csv(file, sep="\t")
        self.data = csv

    def celltypes(self, tissue=None):

        tissues = self.get_tissues()
        if tissue:
            #tissue = tissue.lower()
            assert tissue in tissues, tissue + \
                " does not exist - please pick from: " + " ".join(tissues)
            records = self.data[self.data["organ"] == tissue]
        else:
            records = self.data

        return records["cell type"].unique()

    def marker_matrix(self, celltypes, species=None):
        marker_dict = {}
        records = self.data[self.data["species"]
                            == species] if species else self.data

        for celltype in celltypes:
            celltype_records = records.loc[records["cell type"] == celltype]
            marker_genes = celltype_records["official gene symbol"].unique(
            ).tolist()
            marker_dict[celltype] = marker_genes

        return marker_dict

    def get_tissues(self):
        return [tissue
                for tissue in self.data["organ"].unique() if not isinstance(tissue, float)]


class CellMarker(object):

    def __init__(self):
        pass


panglao = Panglao(file='../downloads/PanglaoDB_markers_17_Dec_2019.tsv')
# all_celltypes = panglao.celltypes()
# print(all_celltypes)

lung_celltypes = panglao.celltypes(tissue="Lungs")
print(lung_celltypes)
matrix = panglao.marker_matrix(lung_celltypes)
print(matrix)


# lung_celltypes2 = panglao.celltypes(tissue="lungs")
# print(lung_celltypes)

# abc_celltypes = panglao.celltypes(tissue="abc")
# print(abc_celltypes)
