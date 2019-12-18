import os
import pandas as pd
import numpy


def get_celltypes(data, celltype_column, tissue_column, tissue=None):
    tissues = get_tissues(data, tissue_column)
    if tissue:
        assert tissue in tissues, tissue + \
            " does not exist - please pick from: " + " ".join(tissues)
        records = data[data[tissue_column] == tissue]
    else:
        records = data

    return records[celltype_column].unique()


def get_tissues(data, column):
    return [tissue
            for tissue in data[column].unique() if not isinstance(tissue, float)]


def get_marker_matrix(data, celltypes, celltype_column, marker_gene_column, species_column, species=None):
    marker_dict = {}
    records = data[data[species_column].isin(
        species)] if species else data

    for celltype in celltypes:
        celltype_records = records.loc[records[celltype_column] == celltype]
        marker_genes = celltype_records[marker_gene_column].unique(
        ).tolist()
        marker_dict[celltype] = marker_genes

    return marker_dict


class Panglao(object):

    def __init__(self, file):
        self.data = pd.read_csv(file, sep="\t")

    def celltypes(self, tissue=None):
        return get_celltypes(self.data, "cell type", "organ", tissue=tissue)

    def marker_matrix(self, celltypes, species=None):
        species = self.get_species_name(species) if species else species
        return get_marker_matrix(self.data, celltypes, "cell type", "official gene symbol", "species", species=species)

    def filter_on_gene_specificity(species):
        pass

    def get_species_name(self, species):
        species = species.lower()
        if species == "human":
            return ["Hs", "Mm Hs"]
        if species == "mouse":
            return ["Mm", "Hs Mm"]


class CellMarker(object):

    def __init__(self, file):
        self.data = pd.read_csv(file, sep="\t")

    def celltypes(self, tissue=None):
        return get_celltypes(self.data, "cellName", "tissueType", tissue=tissue)

    def marker_matrix(self, celltypes, species=None):
        species = self.get_species_name(species) if species else species
        return get_marker_matrix(self.data, celltypes, "cellName", "geneSymbol", "speciesType", species=species)

    def get_species_name(self, species):
        return [species]

    def get_tissues(self):
        return [tissue
                for tissue in self.data["tissueType"].unique() if not isinstance(tissue, float)]


def test_cellMarker():
    cm = CellMarker(file='../downloads/all_cell_markers.txt')

    # all_celltypes = cm.celltypes()
    # print(all_celltypes)

    lung_celltypes = cm.celltypes(tissue="Kidney")
    print(lung_celltypes)
    matrix = cm.marker_matrix(lung_celltypes)
    print(matrix)

    # # lung_celltypes2 = cm.celltypes(tissue="lungs")
    # # print(lung_celltypes)

    # abc_celltypes = cm.celltypes(tissue="abc")
    # print(abc_celltypes)


panglao = Panglao(file='../downloads/PanglaoDB_markers_17_Dec_2019.tsv')
# all_celltypes = panglao.celltypes()
# print(all_celltypes)

test_cellMarker()

lung_celltypes = panglao.celltypes(tissue="Lungs")
print(lung_celltypes)
matrix = panglao.marker_matrix(lung_celltypes)
print(matrix)

# matrix = panglao.marker_matrix(lung_celltypes, species="human")
# print(matrix)


# lung_celltypes2 = panglao.celltypes(tissue="lungs")
# print(lung_celltypes)

# abc_celltypes = panglao.celltypes(tissue="abc")
# print(abc_celltypes)
