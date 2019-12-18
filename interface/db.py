import os
import pandas as pd
import numpy


class Panglao(object):

    def __init__(self, file):
        self.data = pd.read_csv(file, sep="\t")
        

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
        print(self.data["species"].unique())
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

    def filter_on_gene_specificity(species):
        pass

    def get_species_name(species):
        species = species.lower()
        if species == "human":
            pass
        if species == "mouse":
            pass


class CellMarker(object):

    def __init__(self, file):
        self.data = pd.read_csv(file, sep = "\t")


    def celltypes(self, tissue=None):

        tissues = self.get_tissues()
        if tissue:
            #tissue = tissue.lower()
            assert tissue in tissues, tissue + \
                " does not exist - please pick from: " + " ".join(tissues)
            records = self.data[self.data["tissueType"] == tissue]
        else:
            records = self.data

        return records["cellName"].unique()

    def marker_matrix(self, celltypes, species=None):
        marker_dict = {}
        records = self.data[self.data["speciesType"]
                            == species] if species else self.data

        for celltype in celltypes:
            celltype_records = records.loc[records["cellName"] == celltype]
            marker_genes = celltype_records["geneSymbol"].unique(
            ).tolist()
            marker_dict[celltype] = marker_genes

        return marker_dict

    def get_tissues(self):
        return [tissue
                for tissue in self.data["tissueType"].unique() if not isinstance(tissue, float)]




def test_cellMarker():
    cm = CellMarker(file='/Users/abramsd/HACK/cellMarkers.txt')

    all_celltypes = cm.celltypes()
    #print(all_celltypes)

    lung_celltypes = cm.celltypes(tissue="Kidney")
    print(lung_celltypes)
    matrix = cm.marker_matrix(lung_celltypes)
    print(matrix)


    # # lung_celltypes2 = cm.celltypes(tissue="lungs")
    # # print(lung_celltypes)

    # abc_celltypes = cm.celltypes(tissue="abc")
    # print(abc_celltypes)


test_cellMarker()

# panglao = Panglao(file='/Users/abramsd/HACK/cellMarkers.txt')
# # all_celltypes = panglao.celltypes()
# # print(all_celltypes)

# lung_celltypes = panglao.celltypes(tissue="Lungs")
# print(lung_celltypes)
# matrix = panglao.marker_matrix(lung_celltypes)
# print(matrix)


# lung_celltypes2 = panglao.celltypes(tissue="lungs")
# print(lung_celltypes)

# abc_celltypes = panglao.celltypes(tissue="abc")
# print(abc_celltypes)
