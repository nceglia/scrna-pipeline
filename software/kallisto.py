import h5py
import os
import glob
import subprocess
import numpy as np
import time
from scipy.sparse import coo_matrix
from sklearn.preprocessing import normalize
import sys
import collections
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gc
from multiprocessing import Pool
import pickle

from interface.fastqdirectory import FastQDirectory

class Kallisto(object):

    def __init__(self, fastqs, chem="v3"):
        self.output = ".cache/kallisto"
        self.fastqs = fastqs
        self.chem = chem
        self.binary = "kallisto"
        self.index = None
        if not os.path.exists(self.output):
            os.makedirs(self.output)
        self.nthreads = 64
        self.tcc_output = os.path.join(self.output,"tcc")
        if not os.path.exists(self.tcc_output):
            os.makedirs(self.tcc_output)
        self.matrix_ec = os.path.join(self.tcc_output, "matrix.ec")
        self.index = "/reference/Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx"
        self.transcript_to_gene = "/reference/t2g.txt"
        self.bus_output = os.path.join(self.tcc_output,"output.bus")
        self.sorted_bus = os.path.join(self.tcc_output,"sorted.bus")
        self.corrected_bus = os.path.join(self.tcc_output,"corrected.bus")
        self.bus_matrix = os.path.join(self.tcc_output,"matrix.tsv")
        self.bustools = "bustools"
        self.transcripts = os.path.join(self.tcc_output, "transcripts.txt")
        self.tenx_path = os.path.join(self.output, "filtered_feature_bc_matrix")
        if not os.path.exists(self.tenx_path):
            os.makedirs(self.tenx_path)
        self.whitelist = "/reference/10xv3_whitelist.txt"
        self.genes_tsv = os.path.join(self.tenx_path,"genes.genes.txt")
        self.barcodes_tsv = os.path.join(self.tenx_path,"genes.barcodes.txt")
        self.matrix = os.path.join(self.tenx_path,"genes.mtx")

    def run_pseudo(self):
        if not os.path.exists(self.bus_output):
            cmd = [self.binary,"bus","-i",self.index,"-o",self.tcc_output,"-t",str(self.nthreads),"-x","10x{}".format(self.chem)]
            for fastq in self.fastqs.get_fastqs():
                cmd.append(fastq)
            print (" ".join(cmd))
            subprocess.call(cmd)
        assert os.path.exists(self.bus_output)

    def run_bus(self):
        if not os.path.exists(self.corrected_bus):
            print("Running Barcode Correction.")
            cmd = [self.bustools, "correct", "-w", self.whitelist, self.bus_output,"-o",self.corrected_bus]
            subprocess.call(cmd)
        if not os.path.exists(self.sorted_bus):
            print("Sorting Corrected Bus.")
            cmd = [self.bustools, "sort","-T","tmp","-t","64","-m","64G",self.corrected_bus,"-o",self.sorted_bus]
            subprocess.call(cmd)
        if not os.path.exists(self.genes_tsv) or not os.path.exists(self.barcodes_tsv) or not os.path.exists(self.matrix):
            print("Running Transcript Count.")
            cmd = [self.bustools, "count","-o",self.tenx_path + "/genes","-g",self.transcript_to_gene,"-e",self.matrix_ec,"-t",self.transcripts,"--genecounts",self.sorted_bus]
            subprocess.call(cmd)
        assert os.path.exists(self.genes_tsv) and os.path.exists(self.barcodes_tsv) and os.path.exists(self.matrix)

    def transcript_map(self):
        gene_to_transcript = dict()
        genes = open(self.transcript_to_gene,"r").read().splitlines()
        for gene in genes:
            t1, t2, symbol = gene.split()
            gene_to_transcript[t1] = symbol
            gene_to_transcript[t2] = symbol
            gene_to_transcript[t1.rstrip(".1")] = symbol
        return gene_to_transcript

    def count_matrix(self):
        output = open(os.path.join(self.tenx_path,"matrix.mtx"),"w")
        rows = open(self.matrix,"r").read().splitlines()
        for row in rows:
            line = row.split()
            if len(line) == 3:
                new_line = " ".join([line[1],line[0],line[2]])
                output.write(new_line+"\n")
            else:
                output.write(row+"\n")
        output.close()

    def barcodes(self):
        shutil.move(self.barcodes_tsv, os.path.join(self.tenx_path,"barcodes.tsv"))

    def genes():
        output = open(os.path.join(self.tenx_path,"genes.tsv"),"w")
        rows = open(self.genes_tsv,"r").read().splitlines()
        transcripts = transcript_map()
        for row in rows:
            gene = transcripts[row.strip()]
            row = "{}\t{}".format(row.rstrip(".1"), gene)
            output.write(row+"\n")
        output.close()

    def count(self):
        print("Running Kallisto.")
        self.run_pseudo()
        print("Running BUStools.")
        self.run_bus()
        print("Setting up Matrix.")
        self.count_matrix()
        print("Mapping transcripts.")
        self.genes()
        print("Copy valid barcodes.")
        self.barcodes()
        return self.tenx_path


def main():
    sample = "TENX065"
    fastq = "/data/AHT52JBGXB/"
    output = "./"
    fastq_directory = FastQDirectory(fastq, sample, output)

    krunner = Kallisto(fastq_directory)
    tenx_path = krunner.count()
    print(tenx_path)

if __name__ == '__main__':
    main()
