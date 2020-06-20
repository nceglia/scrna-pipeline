import container
import scriptmanager
import martian

__MRO__ = '''
stage BATCH_CORRECTION(
    in  path integrated_seurat,
    in  path image,
    in  string runtime,
    out rds integrated_seurat,
    src py   "stages/batch_correction",
) using (
    threads = 16,
)
'''

def main(args, outs):
    scripts = scriptmanager.ScriptManager()
    script = scripts.batchcorrection()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)
