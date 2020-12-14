import container
import scriptmanager
import martian
import itertools
import os

__MRO__ = '''
stage SUBTYPE_SCORE(
    in string[] celltypes,
    in  map celltype_seurat,
    in  path subtype_yaml,
    in  path image,
    in  string runtime,
    out map subtype_scores,
    src py   "stages/subtype_score",
) split using (
    in  map celltype_seurat,
) using (
    threads = 12,
)
'''

def split(args):
    chunks = []
    subtype_genes = dict()
    subtypes = open("/juno/work/shah/ceglian/chow/transcriptome_analysis/subtypes.csv","r").read().splitlines()
    for subtype in subtypes:
        subtype = subtype.split(",")
        if subtype[0] not in subtype_genes:
            subtype_genes[subtype[0]] = dict()
        subtype_genes[subtype[0]][subtype[1]] = subtype[2].split(";")
    for celltype in args.celltypes:
        if celltype not in subtype_genes: continue
        for subtype, genes in subtype_genes[celltype].items():
            chunk_def = {}
            chunk_def["celltype_obj"] = args.celltype_seurat[celltype]
            chunk_def["subtype"]  = subtype
            chunk_def["celltype"] = celltype
            chunk_def["geneset"] = genes
            chunk_def['__threads'] = 8
            chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    geneset = "{}_{}_genes.txt".format(args.celltype, args.subtype)
    geneset_csv = martian.make_path(geneset)
    output = open(geneset_csv,"w")
    output.write("\n".join(["genes"] + args.geneset))
    output.close()
    args.geneset = geneset_csv
    svg = "{}_{}_subtype.svg".format(args.celltype, args.subtype)
    outs.subtype_scores = martian.make_path(svg)
    scripts = scriptmanager.ScriptManager()
    script = scripts.subtypescores()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    try:
        con.run(script, args, outs)
    except Exception as e:
        return


def join(args, outs, chunk_defs, chunk_outs):
    outs.subtype_scores = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        subtype_name = "{}_{}".format(arg.celltype, arg.subtype)
        if os.path.exists(out.subtype_scores):
            outs.subtype_scores[subtype_name] = {"subtype":subtype_name,"svg":out.subtype_scores}
