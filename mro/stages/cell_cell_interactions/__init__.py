import container
import scriptmanager
import martian
import itertools
import os

__MRO__ = '''
stage CELL_CELL_INTERACTIONS(
    in string[] celltypes,
    in  map celltype_seurat,
    in  path gmt,
    in  path image,
    in  string runtime,
    out map interactions,
    out map network_svg,
    src py   "stages/cell_cell_interactions",
) split using (
    in  map celltype_seurat,
) using (
    threads = 12,
)
'''

def split(args):
    chunks = []
    geneset_oi = dict()
    pathways_of_interest = {"DNA_REPAIR":["Cancer.cell","Cancer.cell"]}
    if not args.gmt or args.gmt == "None":
        return {"chunks": chunks}
    pathways = open(args.gmt,"r").read().splitlines()
    for pathway in pathways:
        pathway = pathway.split()
        geneset_oi[pathway[0].replace("HALLMARK_","")] = pathway[2:]

    for pathway, interactors in pathways_of_interest.items():
        if interactors[0] in args.celltype_seurat and interactors[1] in args.celltype_seurat:
            chunk_def = {}
            chunk_def["sender"]   = interactors[0]
            chunk_def["receiver"] = interactors[1]
            chunk_def["sender_obj"]        = args.celltype_seurat[interactors[0]]
            chunk_def["receiver_obj"]      = args.celltype_seurat[interactors[1]]
            chunk_def["geneset"]       = geneset_oi[pathway.replace("HALLMARK_","")]
            chunk_def["pathway"]       = pathway.replace("HALLMARK_","")
            chunk_def["batch_label"] = ""
            chunk_def['__threads'] = 8
            chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    geneset = "{}_genes.txt".format(args.pathway)
    geneset_csv = martian.make_path(geneset)
    output = open(geneset_csv,"w")
    output.write("\n".join(["genes"] + args.geneset))
    output.close()
    args.geneset = geneset_csv
    png = "{}_{}_{}_interactions.svg".format(args.pathway, args.sender, args.receiver, args.batch_label)
    outs.interactions = martian.make_path(png)
    network_svg = "{}_{}_{}_network.svg".format(args.pathway, args.sender, args.receiver, args.batch_label)
    outs.network_svg = martian.make_path(network_svg)
    weighted_tmp = "{}_{}_{}_weighted_tmp.txt".format(args.pathway, args.sender, args.receiver, args.batch_label)
    outs.weighted_tmp = martian.make_path(weighted_tmp)
    scripts = scriptmanager.ScriptManager()
    script = scripts.cellcellinteractions()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    # try:
    # con.run(script, args, outs)
    # except Exception as e:
    #     pass 

def join(args, outs, chunk_defs, chunk_outs):
    outs.interactions = dict()
    outs.network_svg = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        if os.path.exists(out.interactions):
            if arg.pathway not in outs.interactions:
                outs.interactions[arg.pathway] = dict()
            outs.interactions[arg.pathway][arg.batch_label] = {"sender":arg.sender,"receiver":arg.receiver,"svg":out.interactions}
        if os.path.exists(out.network_svg):
            if arg.pathway not in outs.network_svg:
                outs.network_svg[arg.pathway] = dict()
            outs.network_svg[arg.pathway][arg.batch_label] = {"sender":arg.sender,"receiver":arg.receiver,"svg":out.network_svg}
