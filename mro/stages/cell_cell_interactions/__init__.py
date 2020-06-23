import container
import scriptmanager
import martian
import itertools
import os

__MRO__ = '''
stage CELL_CELL_INTERACTIONS(
    in  string[] celltypes,
    in  map celltype_seurat,
    in  path gmt,
    in  path image,
    in  string runtime,
    out map interactions,
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
    pathways_of_interest = {"HALLMARK_DNA_REPAIR":["Endothelial.cell","Endothelial.cell"],
                            "HALLMARK_INTERFERON_GAMMA_RESPONSE":["B.cell","Endothelial.cell"],
                            "HALLMARK_INFLAMMATORY_RESPONSE":["B.cell","Endothelial.cell"]}
    pathways = open(args.gmt,"r")
    for pathway in pathways:
        pathway = pathway.split()
        geneset_oi[pathway[0].replace("HALLMARK_","")] = pathway[2:]
    for pathway, interactors in pathways_of_interest.items():
        chunk_def = {}
        chunk_def["sender"]   = interactors[0]
        chunk_def["receiver"] = interactors[1]
        chunk_def["sender_obj"]        = args.celltype_seurat[interactors[0]]
        chunk_def["receiver_obj"]      = args.celltype_seurat[interactors[1]]
        chunk_def["geneset"]       = geneset_oi[pathway.replace("HALLMARK_","")]
        chunk_def["pathway"]       = pathway.replace("HALLMARK_","")
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
    png = "{}_{}_{}_interactions.svg".format(args.pathway, args.sender, args.receiver)
    outs.interactions = martian.make_path(png)
    scripts = scriptmanager.ScriptManager()
    script = scripts.cellcellinteractions()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)

def join(args, outs, chunk_defs, chunk_outs):
    outs.interactions = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        if arg.pathway not in outs.interactions:
            outs.interactions[arg.pathway] = []
        if os.path.exists(out.interactions):
            outs.interactions[arg.pathway].append({"sender":arg.sender,"receiver":arg.receiver,"svg":out.interactions})
