import container
import scriptmanager
import martian

__MRO__ = '''
stage CELLTYPE_COMPOSITION(
    in  rds batch_merged_seurat,
    in  path image,
    in  string runtime,
    out png batch_composition,
    out png sample_composition,
    src py   "stages/celltypecomposition",
) using (
    threads = 12,
)
'''

def main(args, outs):
    scripts = scriptmanager.ScriptManager()
    script = scripts.celltypecomposition()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)
