import container
import scriptmanager
import martian

__MRO__ = '''
stage DETECT_DOUBLETS(
    in  map sce,
    in  path image,
    in  string runtime,
    out map csv,
    src py "stages/detect_doublets",
) split using (
    in  map sces,
)
'''

def split(args):
    chunks = []
    for sample in args.sce:
        chunk_def = {}
        chunk_def['sce'] = args.sce[sample]
        chunk_def['sample'] = sample
        chunk_def['__threads'] = 4
        chunk_def['__mem_gb'] = 2
        chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    csv = "{}".format(args.sample)
    outs.csv = martian.make_path(csv)
    scripts = scriptmanager.ScriptManager()
    script = scripts.detectdoublets()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)

def join(args, outs, chunk_defs, chunk_outs):
    outs.csv = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        outs.csv[arg.sample] = out.csv
