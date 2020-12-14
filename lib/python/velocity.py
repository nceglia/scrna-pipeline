import subprocess
import sys
import shutil
import os
import velocyto as vcy
from singlecellexperiment import SingleCellExperiment
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['svg.fonttype'] = 'none'

gtf["HUMAN"] = "/juno/work/shah/reference/transcriptomes/GRCh38/genes/genes.gtf"
mask["HUMAN"] = "/juno/work/shah/ceglian/rnascp/Homo_sapiens.GRCh38.repeat_mask.gtf"
mask["MOUSE"] = "/juno/work/shah/ceglian/chow/transcriptome_analysis/"
gtf["MOUSE"] = "/juno/work/shah/reference/transcriptomes/mm10/genes/genes.gtf"

def run_velocyto(matrix, gtf, mask, directory):
    output_name = directory.split("/")[2]
    outs = "/".join(matrix.split("/")[:-1])
    try:
        os.makedirs(directory)
        os.makedirs(os.path.join(directory,"outs"))
    except Exception as e:
        pass
    bam = os.path.join(outs,"possorted_genome_bam.bam")
    bai = bam + ".bai"
    loomfile = os.path.join(directory,"outs.loom")
    linked_bam = os.path.join(directory,"possorted_genome_bam.bam")
    linked_bai = linked_bam + ".bai"
    filtered_matrix = os.path.join(outs,"filtered_feature_bc_matrix")
    shutil.copytree(filtered_matrix,os.path.join(directory,"outs","filtered_feature_bc_matrix"))
    try:
        os.unlink(linked_bam)
        os.unlink(linked_bai)
    except Exception as e:
        pass
    linked_bam = os.path.join(directory,"outs","possorted_genome_bam.bam")
    linked_bai = linked_bam + ".bai"
    env = os.environ.copy()
    env["PATH"] = "/juno/work/shah/software/samtools-1.10/bin/bin/:" + env["PATH"]
    os.symlink(bam, linked_bam)
    os.symlink(bai, linked_bai)
    cmd = "velocyto run10x -m {mask} {matrix} {gtf} --samtools-threads 16 --samtools-memory 8".format(matrix="/juno"+directory, gtf=gtf, mask=mask)
    subprocess.check_output(cmd.split(),env=env)
    os.unlink(linked_bam)
    os.unlink(linked_bai)
    return os.path.join(outs,"velocyto","{}.loom".format(output_name))

def filter_loom(loomfile, sce, result):
    sce = SingleCellExperiment.fromRData(sce)
    vlm = vcy.VelocytoLoom(loomfile)
    if "Barcode" in sce.colData:
        barcodes = sce.colData["Barcode"]
    else:
        barcodes = sce.colnames
    valid_barcodes = [barcode.split("-")[0] for barcode in barcodes] 
    unfiltered_barcodes = [x.split(":")[1].rstrip("x") for x in vlm.ca["CellID"]]
    filtered_barcodes = [True if barcode in valid_barcodes else False for barcode in unfiltered_barcodes]
    vlm.filter_cells(filtered_barcodes)
    vlm.to_hdf5(result)
    
def analyze_loom(result, fraction_svg):
    vlm = vcy.load_velocyto_hdf5(result)
    vlm.plot_fractions(save2file=fraction_svg)

if __name__=="__main__":
    sample         = sys.argv[1]
    matrix         = sys.argv[2]
    sce         = sys.argv[3]
    loomfile       = sys.argv[4]
    frac_svg      = sys.argv[5]
    directory            = sys.argv[6]
    loom_name      = directory.split("/")[-1]

    #loomfile = run_velocyto(matrix, gtf, mask, directory)
    loomfile = os.path.join("/work/shah/ceglian/rnascp/velocity/","velocity_{}.loom".format(sample))
    result = loomfile.replace(".loom",".hd5")

    filter_loom(loomfile, sce, result)
    analyze_loom(result, frac_svg)
