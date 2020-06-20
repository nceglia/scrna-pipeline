import container
import scriptmanager
import martian

__MRO__ = '''
stage SUBCELLTYPE_MARKERS(
    in  map celltype_seurat
    in  path image,
    in  string runtime,
    out map markers,
    src py   "stages/subcelltype_markers",
) split using (
    in  map celltype_seurat,
) using (
    threads = 12,
)
'''

def split(args):
    chunks = []
    for celltype, path in args.celltype_seurat.items():
        chunk_def = {}
        chunk_def["celltype"] = celltype
        chunk_def['celltype_seurat'] = path
        chunk_def['__threads'] = 8
        chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    csv = "{}_markers.csv".format(args.celltype)
    outs.markers = martian.make_path(csv)
    scripts = scriptmanager.ScriptManager()
    script = scripts.subcelltypemarkers()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)

def join(args, outs, chunk_defs, chunk_outs):
    outs.markers = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        outs.markers[arg.celltype] = out.markers
