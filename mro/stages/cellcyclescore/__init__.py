import container
import scriptmanager
import martian

__MRO__ = '''
stage CELLCYCLE_SCORE(
    in  map qcd_seurat,
    in  path cellcycle_genes,
    in  path image,
    in  string runtime,
    out map qcd_scored_seurat,
    src py   "stages/cellcycle_score",
) split using (
    in  map qcd_seurat,
) using (
    threads = 4,
)
'''

def split(args):
    chunks = []
    for sample in args.qcd_seurat:
        chunk_def = {}
        chunk_def['qcd_seurat'] = args.qcd_seurat[sample]
        chunk_def['sample'] = sample
        chunk_def['__threads'] = 4
        chunk_def['__mem_gb'] = 2
        chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    seurat = "{}_qcd_scored_seurat.rds".format(args.sample)
    outs.qcd_scored_seurat = martian.make_path(seurat)
    sce = "{}_qcd_scored_sce.rds".format(args.sample)
    outs.qcd_scored_sce = martian.make_path(sce)
    scripts = scriptmanager.ScriptManager()
    script = scripts.cellcyclescore()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)

def join(args, outs, chunk_defs, chunk_outs):
    outs.qcd_scored_seurat = dict()
    outs.qcd_scored_sce = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        outs.qcd_scored_seurat[arg.sample] = out.qcd_scored_seurat
        outs.qcd_scored_sce[arg.sample] = out.qcd_scored_sce
