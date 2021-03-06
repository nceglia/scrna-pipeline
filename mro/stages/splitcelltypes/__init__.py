import container
import scriptmanager
import martian

__MRO__ = '''
stage SUBSET_CELLTYPES(
    in  rds batch_merged_seurat,
    in  string[] celltypes,
    in  path image,
    in  string runtime,
    out map celltype_seurat,
    out map cell_umap,
    out map markers,
    src py   "stages/splitcelltypes",
) split using (
    in  string[] celltypes,
) using (
    threads = 12,
)
'''

def split(args):
    chunks = []
    for celltype in args.celltypes:
        chunk_def = {}
        chunk_def['celltype'] = celltype
        chunk_def['__threads'] = 20
        chunk_def['__memgb'] = 4
        chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    rds = "{}_seurat.rds".format(args.celltype)
    outs.celltype_seurat = martian.make_path(rds)
    rds = "{}_sce.rds".format(args.celltype)
    outs.celltype_sce = martian.make_path(rds)
    svg = "{}_celltype.svg".format(args.celltype)
    outs.cell_umap = martian.make_path(svg)
    csv = "{}_clusters.csv".format(args.celltype)
    outs.markers = martian.make_path(csv)
    scripts = scriptmanager.ScriptManager()
    script = scripts.subsetcelltype()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)

def join(args, outs, chunk_defs, chunk_outs):
    outs.celltype_seurat = dict()
    outs.celltype_sce = dict()
    outs.cell_umap = dict()
    outs.markers = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        outs.celltype_seurat[arg.celltype] = out.celltype_seurat
        outs.celltype_sce[arg.celltype] = out.celltype_sce
        outs.cell_umap[arg.celltype] = out.cell_umap
        outs.markers[arg.celltype] = out.markers
