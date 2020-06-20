import container
import scriptmanager
import martian

__MRO__ = '''
stage CELLTYPE_MARKERS(
    in  path yaml,
    in  path image,
    in  string runtime,
    out path marker_csv,
    src py   "stages/celltype_markers",
) using (
    threads = 1,
)
'''

def main(args, outs):
    scripts = scriptmanager.ScriptManager()
    script = scripts.celltypemarkers()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)
