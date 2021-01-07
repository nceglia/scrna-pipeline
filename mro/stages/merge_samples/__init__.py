import container
import scriptmanager
import martian
import scanpy as sc
import pandas
import numpy
import sys

__MRO__ = '''
stage MERGE_SAMPLES(
    in  map scanpy_h5ad,
    in  map batch,
    in  path image,
    in  string runtime,
    out map batch_h5ad,
    src py   "stages/merge_samples",
) split using (
    in  map batch,
) using (
    threads = 4,
    volatile = strict,
)
'''
def split(args):
    chunks = []
    for batch_id, samples in args.batch.items():
        chunk_def = {}
        chunk_def['batched_samples'] = samples
        chunk_def['batch'] = batch_id
        chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    batch_h5ad = "{}.h5ad".format(args.batch)
    merged_tsv = "{}_samples.tsv".format(args.batch)
    outs.batch_h5ad = martian.make_path(batch_h5ad)
    outs.merged_tsv = martian.make_path(merged_tsv)

    adatas = []

    output = open(outs.merged_tsv,"w")
    output.write("sample_id\tpath\n")
    for sample, path in args.scanpy_h5ad.items():
        if sample in args.batched_samples:
            output.write("{}\t{}\n".format(sample,path))
            adata = sc.read(path)
            adata.obs["sample"] = sample
            adatas.append(adata)
    output.close()

    combined = adatas[0].concatenate(*adatas[1:], batch_key="sample")
    sc.pp.filter_cells(combined, min_counts=1)
    sc.pp.calculate_qc_metrics(combined, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pp.normalize_total(combined)
    sc.pp.log1p(combined)
    combined.X = numpy.nan_to_num(combined.X)
    sc.tl.pca(combined, svd_solver='arpack')
    sc.pp.neighbors(combined, n_neighbors=10, n_pcs=50)
    sc.tl.umap(combined)

    combined.write(outs.batch_h5ad.decode())

def join(args, outs, chunk_defs, chunk_outs):
    outs.batch_h5ad = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        outs.batch_h5ad[arg.batch] = out.batch_h5ad
