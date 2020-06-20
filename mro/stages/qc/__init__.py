import container
import scriptmanager
import martian

__MRO__ = '''
stage QC(
    in  path samples,
    in  int  mito,
    in  path image,
    in  string container_runtime,
    out map seurats,
    out map sces,
    src py   "stages/qc",
) split using (
    in  map samples,
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
    seurat = "{}_seurat.rds".format(args.sample)
    sce = "{}_sce.rds".format(args.sample)
    raw_sce = "{}_raw_sce.rds".format(args.sample)
    outs.seurat = martian.make_path(seurat)
    outs.sce    = martian.make_path(sce)
    outs.raw_sce    = martian.make_path(raw_sce)
    scripts = scriptmanager.ScriptManager()
    script = scripts.qc()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)

def join(args, outs, chunk_defs, chunk_outs):
    outs.seurat = dict()
    outs.sce = dict()
    outs.raw_sce = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        outs.seurat[arg.sample] = out.seurat
        outs.sce[arg.sample] = out.sce
        outs.raw_sce[arg.sample] = out.raw_sce
