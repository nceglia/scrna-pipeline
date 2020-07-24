# reference_path: reference_0526/reference.rdata (reference count matrix. same for all patients.)
# gene_path: reference_0526/genes.txt (gene annotation table)
# observation_seurat: patient.rdata (patient level seurat object after CellAssign)
# output_dir (output directory, all infercnv output will be saved in this directory)
# cell_types (use certain cell types for inferCNV, default is Ovarian.cancer.cell)

# example:
# Rscript infercnv_script.R reference_0526/spectrum_infercnv_reference_20200526.rdata reference_0526/genes.txt /work/shah/isabl_data_lake/analyses/19/09/1909/patient.rdata sample_output/

args <- commandArgs(trailingOnly=TRUE)

# reference, count matrix rdata
reference_path <- args[1]
# reference, gene annotation file
gene_path <- args[2]
# seurat object, patient.rdata
observation_seurat <- args[3]
seurat <- readRDS(observation_seurat)
seurat$cell_type <- "Cancer.cell"
saveRDS(seurat, file=observation_seurat)
# directory for all infercnv output
output_dir <- args[4]

# sample cells from input seurat object. The input should be an integer or NA (selecting all cells)
# sample_cell_num <- args[5]
# one or multiple cell_type to be used in infercnv
# if(length(args) >= 5){
#   cell_types <- args[5:length(args)]
# }else{

cell_types <- "Cancer.cell"
# }

print(paste0("reference_path: ", reference_path))
print(paste0("gene_path: ", gene_path))
print(paste0("observation_seurat: ", observation_seurat))
print(paste0("output_dir: ", output_dir))


library(Seurat)

generateInferCNVInput <- function(integratedObj, geneAnnotation, referenceLists, cellTypes = NULL, outputDir, sce = TRUE, sct = FALSE, sampleNumber = NULL){
  geneAnnotation <- geneAnnotation[!duplicated(geneAnnotation$gene_name), ]
  genes <- geneAnnotation$gene_name
  for(i in names(referenceLists)){
    genes <- intersect(genes, rownames(referenceLists[[i]]$expr))
  }
  genes <- intersect(genes, rownames(integratedObj))
  print("here")
  referenceExpr <- NULL
  referenceAnno <- NULL
  for(i in names(referenceLists)){
    referenceExpr <- cbind(referenceExpr, referenceLists[[i]]$expr[genes, ])
    referenceAnno <- rbind(referenceAnno, referenceLists[[i]]$anno)
  }

  if(is.null(cellTypes)){
    observationCells <- colnames(integratedObj)
  }else{
    observationCells <- colnames(integratedObj)[integratedObj$cell_type %in% cellTypes]
  }

  if(!is.null(sampleNumber)){
    if(sampleNumber < length(observationCells)){
      observationCells <- sample(observationCells, sampleNumber)
    }
  }

  if(sce == TRUE){
    observationExpr <- as.matrix(integratedObj@assays$data$counts[genes, observationCells])
  }else{
    if(sct == TRUE){
      observationExpr <- as.matrix(integratedObj@assays$SCT@counts[genes, observationCells])
    }else{
      observationExpr <- as.matrix(integratedObj@assays$RNA@counts[genes, observationCells])
    }

  }
  observationAnno <- data.frame(cell = observationCells, type = "observation")
  geneAnnotation <- geneAnnotation[geneAnnotation$gene_name %in% genes, ]
  geneAnnotation <- geneAnnotation[, c("gene_name", "chromosome", "start", "end")]

  outputMatrix <- as.data.frame(cbind(referenceExpr, observationExpr))
  cellAnnotation <- as.data.frame(rbind(referenceAnno, observationAnno))

  # write
  dir.create(file.path(outputDir), showWarnings = FALSE)
  cellAnnotationPath <- file.path(outputDir, "cellAnnotation.txt")
  geneAnnotationPath <- file.path(outputDir, "geneAnnotation.txt")
  outputMatrixPath <- file.path(outputDir, "countMatrix.txt")
  write.table(cellAnnotation, file = cellAnnotationPath, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  write.table(geneAnnotation, file = geneAnnotationPath, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  write.table(outputMatrix, file = outputMatrixPath, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
}

library(infercnv)
runInferCNV <- function(outputDir, reference = "reference"){
  # mkdir and output the inferCNV input files
  cellAnnotationPath <- file.path(outputDir, "cellAnnotation.txt")
  geneAnnotationPath <- file.path(outputDir, "geneAnnotation.txt")
  outputMatrixPath <- file.path(outputDir, "countMatrix.txt")

  if((!file.exists(cellAnnotationPath)) | (!file.exists(geneAnnotationPath)) | (!file.exists(outputMatrixPath))){
    stop("One or more input files are missing!")
  }

  cellAnnotation <- read.table(cellAnnotationPath, sep = "\t", stringsAsFactors = FALSE)

  reference <- reference[reference %in% cellAnnotation[, 2]]

  infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix= outputMatrixPath,
                                       annotations_file = cellAnnotationPath,
                                       max_cells_per_group=25,
                                       delim = "\t",
                                       gene_order_file = geneAnnotationPath,
                                       ref_group_names = reference)

  infercnv_obj <- infercnv::run(infercnv_obj,
                                cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                min_cells_per_gene = 5,
                                cluster_references = FALSE,
                                out_dir= outputDir,  # dir is auto-created for storing outputs
                                cluster_by_groups=T,   # cluster
                                denoise=T,
                                num_threads=60,
                                output_format = "pdf",
                                HMM=T,
                                HMM_type = "i3",
                                analysis_mode='subclusters')
  #analysis_mode='subclusters',
  #hclust_method='ward.D2',
  #tumor_subcluster_partition_method='random_trees')

  infercnv_obj
}




referenceMatrix <- readRDS(reference_path)
referenceMatrix <- as.matrix(referenceMatrix)


referenceAnnotation <- data.frame(cell = colnames(referenceMatrix), type = "reference")

geneAnnotation <- read.table(gene_path, stringsAsFactors = FALSE)
colnames(geneAnnotation) <- c("gene_name", "gene_id", "chromosome", "start", "end")

seu <- readRDS(observation_seurat)

# generate infercnv input files in output_dir
reference <- list(expr = referenceMatrix, anno = referenceAnnotation)

print("Prepare infercnv input")

generateInferCNVInput(seu, geneAnnotation, list(reference = reference), cell_types, output_dir, sce = FALSE, sct = FALSE)

print("Start inferCNV")

results <- runInferCNV(output_dir)

print("Finish inferCNV")
