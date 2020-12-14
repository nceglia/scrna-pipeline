import container
import scriptmanager
import martian
import os
import shutil

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
    for celltype in args.celltypes.split(","):
        for resolution in (0.1, 0.2, 0.3):
            chunk_def = {}
            chunk_def['celltype'] = celltype
            chunk_def['resolution'] = resolution
            chunk_def['__threads'] = 32
            chunk_def['__memgb'] = 12
            chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    srds = "{}_seurat_{}.rds".format(args.celltype, args.resolution)
    outs.celltype_seurat = martian.make_path(srds)
    rds = "{}_sce_{}.rds".format(args.celltype, args.resolution)
    outs.celltype_sce = martian.make_path(rds)
    svg = "{}_celltype_{}.svg".format(args.celltype, args.resolution)
    outs.cell_umap = martian.make_path(svg)
    csv = "{}_clusters_{}.csv".format(args.celltype, args.resolution)
    outs.markers = martian.make_path(csv)

    markers_tsv = "{}_markers_{}.tsv".format(args.celltype, args.resolution)
    outs.markers_tsv = martian.make_path(markers_tsv)
    cells_tsv = "{}_cells_{}.tsv".format(args.celltype, args.resolution)
    outs.cells_tsv = martian.make_path(cells_tsv)

    scripts = scriptmanager.ScriptManager()
    script = scripts.subsetcelltype()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    try:
        con.run(script, args, outs)
    except Exception as e:
        pass
    try:
        base = outs.markers_tsv.split("RNASCP")[0]
        copied = os.path.join(base, markers_tsv)
        shutil.copyfile(outs.markers_tsv, copied)

        copied = os.path.join(base, cells_tsv)
        shutil.copyfile(outs.cells_tsv, copied)

        copied = os.path.join(base, srds)
        shutil.copyfile(outs.celltype_seurat, copied)
    except Exception as e:
        pass

def join(args, outs, chunk_defs, chunk_outs):
    outs.celltype_seurat = dict()
    outs.celltype_sce = dict()
    outs.cell_umap = dict()
    outs.markers = dict()
    outs.markers_tsv = dict()
    outs.cells_tsv = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        outs.celltype_seurat[arg.celltype] = out.celltype_seurat
        outs.celltype_sce[arg.celltype] = out.celltype_sce
        outs.cell_umap[arg.celltype] = out.cell_umap
        outs.markers[arg.celltype] = out.markers
        outs.markers_tsv[arg.celltype] = out.markers_tsv
        outs.cells_tsv[arg.celltype] = out.cells_tsv