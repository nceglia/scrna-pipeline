# Single Cell RNA-Seq Pipeline #




Pipeline for running single cell rna-seq analysis.

# Running with Docker #

[![Docker Repository on Quay](https://quay.io/repository/nceglia/scrna-pipeline/status "Docker Repository on Quay")](https://quay.io/repository/nceglia/scrna-pipeline)

```
docker run -e "R_HOME=/usr/local/lib/R/" -e "LD_LIBRARY_PATH=/usr/local/lib/R/lib/" -e "PYTHONPATH=$PYTHONPATH:/codebase/SCRNApipeline/" --mount type=bind,source="$(pwd)"/reference,target=/reference --mount type=bind,source="$(pwd)"/results,target=/results --mount type=bind,source="$(pwd)"/data,target=/data --mount type=bind,source="$(pwd)/runs",target=/runs -w="/runs" -t nceglia/scrna-pipeline:v1.0.0 run_vm --sampleid test --build test
```

# Running with Singularity #

[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/2919)

# Setup #

1. Required arguments.
- sampleid (test)
- build (GRCh38,mm10,test)

Create a job directory with the directory structure including jobs/ data/ reference/.
```
mkdir somedirectory/test-run
mkdir somedirectory/test-run/results
mkdir somedirectory/test-run/data
mkdir somedirectory/test-run/reference
mkdir somedirectory/test-run/runs
```

Run docker on VM (not lsf).


Main results are stored in mounted volumes and loaded into azure.
1. CellRanger Results in `cellrangerv3`
2. Aligned Bams in `bams`
2. QC'd SCEs in `rdatav3`
3. HTML report, figures, cellassign, clonealign, scvis compressed in `results`


# SCRNA Viz json

## Overview
 - Found in **scrnaviz** container in **scrnadata** storage account.
Right now, just an example. I will figure out how to index patient ID for subsequent runs.

- Download example: *patient_data.tar.gz*

-- patient_data.json

```
import json
data = json.loads(open("patient_data.json","r").read())
```

-- *htmls*
Just a directory of cellranger web summary htmls.

### Fields
```
fields = data.keys()
```

1. Sample Names (Sites)
 - Samples should start with "*Sample*" right now.
```
samples = list(filter(lambda x: "Sample" in x, fields))
right_ovary = data["SampleSite"]
right_ovary_cells = right_ovary["celldata"]["Barcode"] # All barcodes for sample
right_ovary_genes = right_ovary["genedata"]["Symbol"] # All genes for sample
```

2. Statistics
 - These values were taken from the web summary or qc R scripts.  They are currently strings :(.
They are all fairly important, try to compare all values across sites.
```
stats_by_samples = data["statistics"]
right_ovary_stats = stats_by_samples["SampleSite"]
cellranger_kit = right_ovary_stats["Chemistry"]
estimated_number_of_cells = int(right_ovary_stats["Estimated Number of Cells"])
cells_retained_at_default_value_for_qc = int(right_ovary_stats["Mito10"])
assert estimated_number_of_cells > cells_retained_at_default_value_for_qc, "We keep less cells than cellranger finds."
saturation = right_ovary_stats["Sequencing Saturation"]
```

3. Cellassign Rho
- This includes all the expected cell types and the marker genes associated.
- Will probably want to combine celltype with marker gene expression (heat).
```
rho = data["rho"]
cell_types = list(rho.keys())
all_marker_genes = list(rho.values())
plasma_cell_marker_genes = rho["Cytotoxic T cells"]
```

### Sample (Site) specific data

1. Matrix
- Count expression (heat) matrix in cells by genes.
- Matrix is in a sparse format, if the cell does not have expression for a specific gene, will raise KeyError.

```
matrix = right_ovary["matrix"]
count = matrix["AAACCCAAGGCATGCA-1"]["CDK11A"]
try:
  count = matrix["AAACCCAAGGCATGCA-1"]["MIR1302-10"] #No expression
except KeyError:
  count = 0
```

2. tSNE
- Coordinates for QC'd cells (Mito10).
```
tsne = right_ovary["tsne"]
x_coordinate, y_coordinate = tsne["AAACCCAAGGCATGCA-1"]
```

3. Cell Types
- Cell types for QC'd cells (Mito10).
```
celltype = right_ovary["cellassign"]["AAACCCAAGGCATGCA-1"]
assert celltype in celltypes, "This celltype should be expected in the rho matrix."
# Get the marker gene expression in this cell for the given cell type.
marker_genes = rho[celltype]
counts_by_gene = matrix["AAACCCAAGGCATGCA-1"]
expression_subset = dict((gene, counts_by_gene[gene]) for gene in marker_genes if gene in counts_by_gene)
```

4. Clusters
- Cluster assignments for QC'd cells from Leiden algorithm.
```
cluster = right_ovary["clusters"]["AAACCCAAGGCATGCA-1"]
```

5. Cell Data

```
celldata = right_ovary["celldata"]
total_counts_to_be_binned = celldata["total_counts"] #make me a histo or a violin and compare across sites :)
total_counts_to_be_binned = celldata["total_features_by_counts"]
total_counts_mito_to_be_binned = celldata["total_counts_mito"]
pct_counts_mito_to_be_binned = celldata["pct_counts_mito"] #Use me to show finer resolution on whether 10% is a good QC :)
total_counts_ribo_to_be_binned = celldata["total_counts_ribo"]
```


##### Things to add:
1) Umap and scvis
2) ... will think more ...

##### Things to fix:
1) Create issues as needed to fix the format of the data.
