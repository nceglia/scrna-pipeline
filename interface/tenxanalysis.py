import os
import collections
import glob
import scanpy.api as sc
from sklearn import cluster
import pandas
import subprocess
from software.tenx import TenX
from utils.config import Configuration
from utils.storage import TenxDataStorage
import matplotlib.pyplot as plt
import numpy
import shutil
import gzip
import tarfile
import pyparsing as pp
import pickle
from interface.singlecellexperiment import SingleCellExperiment
from scipy import io, sparse


config = Configuration()

class TenxAnalysis(object):

    def __init__(self, directory):
        self.loaded = False
        self.directory = directory
        self.path = directory

    def load(self):
        if not os.path.exists(self.directory):
            self.sample = self.directory
            if not os.path.exists(".cache/{}".format(self.sample)):
                cloud_storage = TenxDataStorage(self.sample)
                self.directory = cloud_storage.download()
            else:
                self.directory = ".cache/{}".format(self.sample)
        self.path = self.directory
        v3_path_raw = self.raw_gene_bc_matrices = os.path.join(self.path, 'raw_feature_bc_matrix')
        v2_path_raw = self.raw_gene_bc_matrices = os.path.join(self.path, 'raw_gene_bc_matrices')
        if os.path.exists(v3_path_raw):
            self.raw_gene_bc_matrices = v3_path_raw
            self.detected_version = "v3"
        elif os.path.exists(v2_path_raw):
            self.raw_gene_bc_matrices = v2_path_raw
            self.detected_version = "v2"
        elif os.path.exists(v3_path_raw + "_mex"):
            self.raw_gene_bc_matrices = v3_path_raw + "_mex"
            self.detected_version = "v3"
        elif os.path.exists(v2_path_raw + "_mex"):
            self.raw_gene_bc_matrices = v2_path_raw + "_mex"
            self.detected_version = "v2"
        else:
            print("No Raw Matrices folder found -- Check dir name (raw_feature_bc_matrix or raw_gene_bc_matrices)")
        v3_path_filtered = os.path.join(self.path, 'filtered_feature_bc_matrix')
        v2_path_filtered = os.path.join(self.path, 'filtered_gene_bc_matrices')
        if os.path.exists(v3_path_filtered):
            self.filtered_gene_bc_matrices = v3_path_filtered
            self.detected_version = "v3"
        elif os.path.exists(v2_path_filtered):
            self.filtered_gene_bc_matrices = v2_path_filtered
            self.detected_version = "v2"
        elif os.path.exists(v3_path_filtered + "_mex"):
            self.filtered_gene_bc_matrices = v3_path_filtered + "_mex"
            self.detected_version = "v3"
        elif os.path.exists(v2_path_filtered + "_mex"):
            self.filtered_gene_bc_matrices = v2_path_filtered + "_mex"
            self.detected_version = "v2"
        else:
            print("No Filtered Matrices folder found -- Check dir name (filtered_feature_bc_matrix or filtered_gene_bc_matrices)")
        self.clustering = os.path.join(self.path, 'analysis/clustering')
        self.matrix = os.path.join(self.path, "")
        self.projection = os.path.join(self.path, 'analysis/pca/10_components/projection.csv')
        self.cellranger_tsne = os.path.join(self.path, 'analysis/tsne/2_components/projection.csv')
        self.summary = os.path.join(self.path,"web_summary.html")
        self.metrics_summary = os.path.join(self.path, "metrics_summary.csv")
        self.top_level = "/".join(self.path.split("/")[:-3])

        self.baseobj = "sce.rdata"
        self.qcdobj = "qcdsce.rdata"
        self.rdata = os.path.join(self.directory,self.baseobj)
        self.qcdrdata = os.path.join(self.directory, self.qcdobj)

    def decompress(self,gzipped,extracted):
        with gzip.open(gzipped, 'rb') as f_in:
            with open(extracted, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    def extract(self):
        def check_and_decompress(gzipped,flat):
            if os.path.exists(gzipped) and not os.path.exists(flat):
                self.decompress(gzipped, flat)
        try:
            filtered = self.filtered_matrices()
            self.gzipped_filtered_barcodes = os.path.join(filtered, "barcodes.tsv.gz")
            self._filtered_barcodes = self.gzipped_filtered_barcodes.replace(".gz","")
            check_and_decompress(self.gzipped_filtered_barcodes,self._filtered_barcodes)
            self.gzipped_filtered_matrices = os.path.join(filtered, "matrix.mtx.gz")
            self._filtered_matrices = self.gzipped_filtered_matrices.replace(".gz","")
            check_and_decompress(self.gzipped_filtered_matrices,self._filtered_matrices)
            self.gzipped_filtered_genes = os.path.join(filtered, "features.tsv.gz")
            self.filtered_genes = self.gzipped_filtered_genes.replace("features","genes").replace(".gz","")
            check_and_decompress(self.gzipped_filtered_genes,self.filtered_genes)
        except Exception as e:
            pass
        try:
            raw = self.raw_matrices()
            self.gzipped_raw_barcodes = os.path.join(raw, "barcodes.tsv.gz")
            self.raw_barcodes = self.gzipped_raw_barcodes.replace(".gz","")
            check_and_decompress(self.gzipped_raw_barcodes,self.raw_barcodes)
            self.gzipped_raw_matrices = os.path.join(raw, "matrix.mtx.gz")
            self._raw_matrices = self.gzipped_raw_matrices.replace(".gz","")
            check_and_decompress(self.gzipped_raw_matrices,self._raw_matrices)
            self.gzipped_raw_genes = os.path.join(raw, "features.tsv.gz")
            self.raw_genes = self.gzipped_raw_genes.replace("features","genes").replace(".gz","")
            check_and_decompress(self.gzipped_raw_genes,self.raw_genes)
        except Exception as e:
            pass

    @property
    def chemistry(self):
        return self._chemistry

    @chemistry.getter
    def chemistry(self):
        rows = open(self.summary,"r").read().splitlines()
        for i, row in enumerate(rows):
            if "Chemistry" in row:
                break
        chem = rows[i+1].strip().replace("<td>","").replace("</td>","")
        return chem

    @property
    def metrics(self):
        return self._metrics

    def set_integrated(self, integrated):
        integrated_pkl = os.path.join(self.directory, "integrated.pkl")
        pickle.dump(integrated,open(integrated_pkl,"wb"))

    def get_integrated(self):
        integrated_pkl = os.path.join(self.directory, "integrated.pkl")
        return pickle.load(open(integrated_pkl,"rb"))

    def set_corrected(self, corrected):
        corrected_pkl = os.path.join(self.directory, "corrected.pkl")
        pickle.dump(corrected,open(corrected_pkl,"wb"))

    def get_corrected(self, corrected):
        corrected_pkl = os.path.join(self.directory, "corrected.pkl")
        return pickle.load(open(corrected_pkl,"rb"))

    @metrics.getter
    def metrics(self):
        rows = open(self.metrics_summary,"r").read().splitlines()
        header = rows.pop(0)
        header = pp.commaSeparatedList.parseString(header).asList()
        stats = rows.pop(0)
        stats = pp.commaSeparatedList.parseString(stats).asList()
        assert len(header) == len(stats), "{} - {}".format(len(header),len(stats))
        return dict(zip(header,stats))

    def bus_finalize(self):
        base = "/".join(self.path.split("/")[:-2])
        sample = self.path.split("/")[-2] + ".tar.gz"
        self.outstarball = os.path.join(base, sample)
        with tarfile.open(self.outstarball, "w:gz") as tar:
            tar.add(self.path, arcname=os.path.basename(self.path))

    def finalize(self):
        base = "/".join(self.path.split("/")[:-2])
        print("Outs: ", base)
        bamdir = os.path.join(base,"bams")
        try:
            os.makedirs(bamdir)
        except Exception as e:
            pass
        bams = glob.glob(self.path + "/pos*")
        for bam in bams:
            shutil.move(bam,bamdir)
        self.bamtarball = os.path.join(base, "bams.tar.gz")
        sample = self.path.split("/")[-3] + ".tar.gz"
        self.outstarball = os.path.join(base, sample)
        with tarfile.open(self.outstarball, "w:gz") as tar:
            tar.add(self.path, arcname=os.path.basename(self.path))
        with tarfile.open(self.bamtarball, "w:gz") as tar:
            tar.add(bamdir, arcname=os.path.basename(bamdir))

    def filtered_sce(self):
        return TenX.read10xCountsFiltered(tenx,rdata)

    @staticmethod
    def make_10x_output(adata, path):
        if not os.path.exists(path):
            os.makedirs(path)

        output = open(os.path.join(path, "barcodes.tsv"),"w")
        for barcode in adata.obs.index:
            output.write(barcode+"\n")
        output.close()

        output = open(os.path.join(path,"genes.tsv"),"w")
        for gene in adata.var.index:
            output.write(gene+"\t"+gene+"\n")
        output.close()

        output = str(os.path.join(path, "matrix.mtx"))
        io.mmwrite(output,adata.X.T)


    def adata(self, scepath, subset=None):
        if scepath is None:
            scepath = self.rdata
        scepath = os.path.abspath(scepath)
        print(scepath)
        sce = SingleCellExperiment.fromRData(scepath)
        return self.create_scanpy_adata(sce, subset=subset)

    def sce(self):
        return SingleCellExperiment.fromRData(self.rdata)

    def bam_tarball(self):
        return self.bamtarball

    def outs_tarball(self):
        return self.outstarball

    def filtered_matrices(self):
        if self.detected_version == "v3":
            _base = self.filtered_gene_bc_matrices
            _features = os.path.join(_base, "features.tsv")
            _genes = os.path.join(_base,"genes.tsv")
            if os.path.exists(_features) and not os.path.exists(_genes):
                shutil.copyfile(_features,_genes)
            return self.filtered_gene_bc_matrices
        else:
            build_filtered = os.path.join(self.filtered_gene_bc_matrices, config.build + "/")
            if os.path.exists(build_filtered):
                return build_filtered
            else:
                possible = glob.glob(os.path.join(self.filtered_gene_bc_matrices, "*"))
                if len(possible) == 1:
                    return possible[0]
        raise ValueError("No filtered matrices found")

    def raw_matrices(self):
        if self.detected_version == "v3":
            return self.raw_gene_bc_matrices
        else:
            build_raw = os.path.join(self.raw_gene_bc_matrices, config.build + "/")
            if os.path.exists(build_raw):
                return build_raw
            else:
                possible = glob.glob(os.path.join(self.raw_gene_bc_matrices, "*"))
                if len(possible) == 1:
                    return possible[0]
        raise ValueError("No filtered matrices found")

    def filtered_barcodes(self):
        barcode_file = os.path.join(self.filtered_matrices(),"barcodes.tsv")
        return open(barcode_file,"r").read().splitlines()

    def raw_barcodes(self):
        barcode_file = os.path.join(self.raw_matrices(),"barcodes.tsv")
        return open(barcode_file,"r").read().splitlines()

    def raw_genes(self):
        genes_file = os.path.join(self.raw_matrices(),"genes.tsv")
        return dict([row.split() for row in open(genes_file,"r").read().splitlines()])

    def filtered_genes(self, as_list=False):
        genes_file = os.path.join(self.filtered_matrices(),"genes.tsv")
        if as_list:
            return [row.split()[0] for row in open(genes_file,"r").read().splitlines()]

        return dict([row.split() for row in open(genes_file,"r").read().splitlines()])

    def summary(self):
        web_summary = os.path.join(self.path, "web_summary.html")
        return web_summary

    def molecules_h5(self):
        return os.path.join(self.path,"molecule_info.h5")


    def get_genes(self, sce):
        _transcripts = sce.rowData["hgnc_symbol"]
        try:
            adata = sc.read_10x_h5(self.filtered_h5(), genome=config.build)
        except Exception:
            adata = sc.read_10x_mtx(self.filtered_matrices())
        transcripts = []
        for symbol in transcripts:
            if symbol not in adata.var.index:
                symbol = symbol.replace(".","-")
                if symbol not in adata.var.index:
                    symbol = symbol.split("-")
                    symbol = "-".join(symbol[:-1]) + ".{}".format(symbol[-1])
                    if symbol not in adata.var.index:
                        symbol = symbol.split(".")[0]
            transcripts.append(symbol)
        return _transcripts

    def gene_map(self, sce, original=False):
        _transcripts = sce.rowData["Symbol"]
        try:
            adata = sc.read_10x_h5(self.filtered_h5(), genome=config.build)
        except Exception:
            adata = sc.read_10x_mtx(self.filtered_matrices())
        transcripts = {}
        for symbol in _transcripts:
            original = symbol
            if symbol not in adata.var.index:
                symbol = symbol.replace(".","-")
                if symbol not in adata.var.index:
                    symbol = symbol.split("-")
                    symbol = "-".join(symbol[:-1]) + ".{}".format(symbol[-1])
                    if symbol not in adata.var.index:
                        symbol = symbol.split(".")[0]
            if original:
                transcripts[original] = symbol
            else:
                transcripts[symbol] = original
        return transcripts

    def create_scanpy_adata_basic(self, assay="counts", sample_key=None):
        adata = sc.read_10x_mtx(self.filtered_matrices(), make_unique=True)
        return adata


    def create_scanpy_adata(self, sce, fast_load=True, assay="counts", high_var = False, subset=None):
        barcodes = sce.colData["Barcode"]
        _transcripts = sce.rowData["hgnc_symbol"]
        adata = sc.read_10x_mtx(self.filtered_matrices(), make_unique=True)
        adata.var_names_make_unique()
        adata.obs_names_make_unique()
        adata.barcodes = pandas.read_csv(os.path.join(self.filtered_matrices(),'barcodes.tsv'), header=None)[0]
        adata = adata[barcodes,:]
        return adata

    def get_scvis_dimensions(self, embedding_file):
        if embedding_file is None:
            raise AssertionError("scvis requires embedding file.")
        rows = open(embedding_file,"r").read().splitlines()
        header = rows.pop(0)
        embedding = []
        for row in rows:
            row = list(map(float, row.split("\t")[1:]))
            embedding.append(row)
        return numpy.array(embedding).reshape(2,len(rows))

    def clusters(self, sce, pcs=50, copy = False):
        adata = self.create_scanpy_adata_basic()
        adata = sc.tl.pca(adata, copy=True)
        adata = sc.pp.neighbors(adata, copy=True)
        sc.tl.leiden(adata)
        if not copy:
            return adata.obs["leiden"]
        else:
            return adata

    def markers(self, sce, n=100, pcs=50):
        adata = self.clusters(sce, pcs=pcs,copy=True)
        sc.tl.rank_genes_groups(adata, 'leiden', n_genes=n)
        markers = []
        for genes, pvals in zip(adata.uns["rank_genes_groups"]["names"],adata.uns["rank_genes_groups"]["pvals_adj"]):
            for gene, pval in zip(genes,pvals):
                if pval < 0.01:
                    markers.append(gene)
        return list(set(markers))

    def markers_by_clusters(self, sce, rep=None, pcs=50, embedding_file=None):
        adata = self.clusters(sce, pcs=pcs, embedding_file=embedding_file, copy=True)
        sc.tl.rank_genes_groups(adata, 'leiden')
        sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, save='_{}_{}.png'.format(rep,pcs))
        return adata.uns

    def filtered_h5(self):
        return os.path.join(self.path, "filtered_gene_bc_matrices_h5.h5")

    def raw_h5(self):
        return os.path.join(self.path, "raw_gene_bc_matrices_h5.h5")

    def filtered_mtx(self, genes, barcodes):
        matrix = os.path.join(self.filtered_matrices(),"matrix.mtx")
        rows = open(matrix,"r").read().splitlines()
        print(rows.pop(0))
        print(rows.pop(0))
        lgenes, lcells, lnzero = rows.pop(0).split()
        print(lgenes, lcells)
        sparse_matrix = collections.defaultdict(dict)
        count = 0
        for row in rows:
            row, col, val = row.split()
            gene = genes[int(row)-1]
            barcode = barcodes[int(col)-1]
            sparse_matrix[gene][barcode] = int(val)
            count += 1
        return sparse_matrix

    @staticmethod
    def read_mtx_csv(path):
        data = dict()
        genes = open(path,"r").read()
        barcodes = genes.pop(0).split(",")[1:]
        print(len(barcodes))
        for gene in genes:
            umi_counts = list(map(int, gene[1:]))
            data[gene[0]] = dict(zip(barcodes,umi_counts))
        return data

    def __add__(self, other):

        self_barcodes = self.filtered_barcodes()
        other_barcodes = other.filtered_barcodes()

        self_genes = self.filtered_genes(as_list=True)
        other_genes = self.filtered_genes(as_list=True)
        self_gene_dict = self.filtered_genes()
        other_gene_dict = other.filtered_genes()

        self_sparse_mtx = self.filtered_mtx(self_genes, self_barcodes)
        other_sparse_mtx = self.filtered_mtx(other_genes, other_barcodes)

        barcodes = set(self_barcodes).union(set(other_barcodes))
        print("Reading Mtx")
        union_genes = set(self_sparse_mtx.keys()).union(set(other_sparse_mtx.keys()))
        matrix = dict()
        nonzero = 0
        print ("Processing Genes")
        from tqdm import tqdm
        for gene in tqdm(union_genes):
            self_umis = self_sparse_mtx[gene]
            other_umis = other_sparse_mtx[gene]
            row = dict()
            for barcode in barcodes:
                self_count = 0
                other_count = 0
                if barcode in self_umis: self_count = self_umis[barcode]
                if barcode in other_umis: other_count = other_umis[barcode]
                if self_count + other_count == 0: continue
                row[barcode] = self_count + other_count
                nonzero += 1
            matrix[gene] = row
        combined_mtx = os.path.join(self.path, "combined_gene_bc_matrices")
        try:
            os.makedirs(combined_mtx)
        except OSError:
            pass
        mtx = os.path.join(combined_mtx, "matrix.mtx")
        genes = os.path.join(combined_mtx, "genes.tsv")
        cells = os.path.join(combined_mtx, "barcodes.tsv")
        output = open(mtx,"w")
        output.write("%%MatrixMarket matrix coordinate integer general\n%\n")
        output.write("{} {} {}\n".format(len(union_genes),len(barcodes),nonzero))
        all_genes = matrix.keys()
        row = 0
        col = 0
        use_genes = []
        use_barcodes = []
        for gene in tqdm(union_genes):
            for barcode in barcodes:
                try:
                    output.write("{} {} {}\n".format(row, col, matrix[gene][barcode]))
                    row += 1
                    col += 1
                except Exception as e:
                    continue
        output.close()
        output = open(genes,"w")
        for gene in use_genes:
            symbol = self_gene_dict.get(gene,"---")
            if symbol == "---": symbol = other_gene_dict.get(gene,"---")
            output.write("{} {}\n".format(gene,symbol))
        output.close()
        output = open(cells, "w")
        for barcode in use_barcodes:
            output.write(barcode +"\n")
        output.close()
        return combined_mtx
