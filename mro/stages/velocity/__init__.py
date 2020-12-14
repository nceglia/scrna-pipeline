import container
import scriptmanager
import martian
import os

__MRO__ = '''
stage VELOCITY(
    in  path bam_inventory,
    in  map samples,
    in  map raw_sce,
    in  path image,
    in  string runtime,
    out map frac_svg,
    src py "stages/velocity",
) split using (
    in  map samples,
) using (
    threads = 4,
)
'''

def split(args):
    chunks = []
    outs_dirs = dict()
    if not args.bam_inventory or args.bam_inventory == "None":
        return {"chunks": chunks}
    rows = open(args.bam_inventory,"r").read().splitlines()
    header = rows.pop(0)
    for row in rows:
        sample, outs = row.split(",")
        outs_dirs[sample.split("-")[0]] = outs
    for sample in args.samples:
        if sample not in outs_dirs: continue
        chunk_def = {}
        chunk_def['sample'] = sample
        chunk_def['matrix'] = outs_dirs[sample.split("-")[0]]
        chunk_def['seurat'] = args.samples[sample]
        chunk_def['sce']    = args.sce[sample]
        chunk_def['__threads'] = 16
        chunk_def['__mem_gb'] = 8
        chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    frac_svg = "{}_fractions.svg".format(args.sample)
    outs.frac_svg = martian.make_path(frac_svg)
    loomfile = "{}.loom".format(args.sample)
    outs.looms = martian.make_path(loomfile)
    outputdir = "velocity_{}".format(args.sample)
    outs.directory = martian.make_path(outputdir)
    scripts = scriptmanager.ScriptManager()
    script = scripts.velocity()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)

def join(args, outs, chunk_defs, chunk_outs):
    outs.frac_svg = dict()
    outs.looms    = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        if os.path.exists(out.looms):
            outs.looms[arg.sample] = out.looms
        if os.path.exists(out.frac_svg):
            outs.frac_svg[arg.sample] = out.frac_svg
        