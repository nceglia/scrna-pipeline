# Single Cell RNA-Seq Pipeline #


Pipeline for running single cell rna-seq analysis.

![workflow](resources/rnascp.png)


## Setup ##

 - Download Martian: https://martian-lang.org/getting-started/.
 - Make sure Martian bin directory is in path.
 - source source.sh

## Inputs ##

#### Inventory ####

CSV file with paths to market matrix files produced by tenx cellranger, seurat, etc.

| SampleID    | BatchID     | Path
| ----------- | ----------- | ---------------- |
| Sample1     | Patient1    | "path/to/matrix"
| Sample2     | Patient2    | "path/to/matrix"

#### GMT ####

Geneset GMT file for pathway analysis (see resources/h.all.v7.0.symbols.gmt)

#### Cell type Marker YAML ####

Cell Type Marker Matrix YAML (see resources/hgsc_v5_major.yaml).

```
T.Cell:
  - Gene1
  - Gene2
  - Gene3
 B.Cell:
  - Gene1
  - Gene3
```

#### Subtype CSV ####

CSV in the format (genes separated by semicolon).

| Subtype     | CellType    | Marker Genes
| ----------- | ----------- | ---------------- |
| T.Cell      | Naive       | Gene4;Gene5;Gene6
| T.Cell      | Memory      | Gene6;Gene10


## Running ##

```
rnascp --project awesome --mito 25 --doublet 0.25 --image "nceglia/scrna-pipeline:latest" --runtime singularity --mode lsf --gmt resources/h.all.v7.0.symbols.gmt --yaml resources/hgsc_v5_major.yaml --jobs 40 --mempercore=12
```




