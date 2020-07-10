import container
import scriptmanager
import martian
import itertools
import os

__MRO__ = '''
stage CELLTYPE_VELOCITY(
    in  string[] celltypes,
    in  map celltype_sce,
    in  map looms,
    in  path image,
    in  string runtime,
    out map velocity,
    src py   "stages/celltype_velocity",
) split using (
    in  map celltype_sce,
) using (
    threads = 12,
)
'''


def split(args):
    chunks = []
    for celltype in args.celltypes:
        chunk_def = {}
        chunk_def["sce"] = args.celltype_sce[celltype]
        chunk_def["celltype"] = celltype
        chunk_def["looms"] = args.looms
        chunk_def['__threads'] = 8
        chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    velocity_plot = "{}_velocity.svg".format(args.celltype)
    velocity_plot = martian.make_path(velocity_plot)
    outs.velocity = velocity_plot

    loom = "{}.loom".format(args.celltype)
    loom = martian.make_path(loom)
    outs.celltype_loom = loom
    args.reduction = "celltypeumap"

    scripts = scriptmanager.ScriptManager()
    script = scripts.celltypevelocity()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)

def join(args, outs, chunk_defs, chunk_outs):
    outs.velocity = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        outs.velocity[arg.celltype] = {"svg": out.velocity}
        
