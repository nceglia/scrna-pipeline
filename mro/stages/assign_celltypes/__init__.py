import container
import scriptmanager
import martian
import os

__MRO__ = '''
stage ASSIGN_CELLTYPES(
    in  map batch_h5ad,
    in  path marker_csv,
    in  path image,
    in  string runtime,
    out map celltype_csv,
    src py   "stages/assign_celltypes",
) split using (
    in  map batch_h5ad,
) using (
    threads = 30,
    volatile = strict,
)
'''
def split(args):
    chunks = []
    for batch_id, path in args.merged_seurat.items():
        chunk_def = {}
        chunk_def['merged_seurat'] = path
        chunk_def['batch'] = batch_id
        chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    celltypes_csv  = "{}_celltypes.csv".format(args.batch)
    outs.celltypes_csv    = martian.make_path(celltypes_csv)
    scripts = scriptmanager.ScriptManager()
    script = scripts.assigncelltypes()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)

def join(args, outs, chunk_defs, chunk_outs):
    outs.batch_csv = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        outs.batch_csv[arg.batch] = out.batch_csv
