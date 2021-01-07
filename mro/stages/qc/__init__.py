import container
import scriptmanager
import martian
import scanpy as sc
import scrublet

__MRO__ = '''
stage RUN_QC(
    in  map samples,
    in  int mito,
    in  int ncounts,
    in  int nfeatures,
    in  int mincells,
    in  float scrublet_score,
    in  string cell_cycle_genes,
    in  path image,
    in  string runtime,
    out map scanpy_h5ad,
    src py "stages/qc",
) split using (
    in  map samples,
) using (
    threads = 1,
    volatile = strict,
)
'''

def split(args):
    chunks = []
    for sample in args.samples:
        chunk_def = {}
        chunk_def['sample'] = sample
        chunk_def['matrix'] = args.samples[sample]
        chunk_def['__threads'] = 4
        chunk_def['__mem_gb'] = 2
        chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    scanpy_h5ad = "{}.h5ad".format(args.sample)
    outs.scanpy_h5ad = martian.make_path(scanpy_h5ad)
    adata = sc.read_10x_mtx(args.matrix, var_names='gene_symbols', cache=True)

    scrub = scrublet.Scrublet(adata.X)
    doublet_scores, _ = scrub.scrub_doublets(n_prin_comps=10)

    adata.obs["doublet_score"] = doublet_scores
    adata = adata[adata.obs["doublet_score"] < args.scrublet_score]

    sc.pp.filter_cells(adata, min_genes=args.nfeatures)
    sc.pp.filter_genes(adata, min_cells=args.mincells)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs.total_counts > args.ncounts, :]
    adata = adata[adata.obs.pct_counts_mt < args.mito, :]
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
    sc.tl.umap(adata)

    cell_cycle_genes = [x.strip() for x in open(args.cell_cycle_genes)]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    sc.pp.regress_out(adata, ['total_counts','pct_counts_mt','S_score','G2M_score'])
    sc.pp.scale(adata)

    adata.write(outs.scanpy_h5ad.decode())

def join(args, outs, chunk_defs, chunk_outs):
    outs.scanpy_h5ad = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        outs.scanpy_h5ad[arg.sample] = out.scanpy_h5ad
