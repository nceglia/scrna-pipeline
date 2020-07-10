import container
import scriptmanager
import martian

__MRO__ = '''
stage FILTER_DOUBLETS(
    in  map seurat,
    in  map csv,
    in  float score,
    in  path image,
    in  string runtime,
    out map qcd_seurat,
    src py   "stages/filter_doublets",
) split using (
    in  map seurat,
    in  map csv,
)
'''
def split(args):
    chunks = []
    for sample in args.seurat:
        chunk_def = {}
        chunk_def['sample'] = sample
        chunk_def['seurat'] = args.seurat[sample]
        chunk_def['csv'] = args.csv[sample]
        chunk_def['__threads'] = 4
        chunk_def['__mem_gb'] = 2
        chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    qcd_seurat = "{}_qcd_seurat.rds".format(args.sample)
    outs.qcd_seurat = martian.make_path(qcd_seurat)
    scripts = scriptmanager.ScriptManager()
    script = scripts.filterdoublets()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)

def join(args, outs, chunk_defs, chunk_outs):
    outs.qcd_seurat = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        outs.qcd_seurat[arg.sample] = out.qcd_seurat
