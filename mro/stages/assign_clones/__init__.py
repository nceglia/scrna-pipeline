import container
import scriptmanager
import martian
import os

__MRO__ = '''
stage ASSIGN_CELLTYPES(
    in  map annotated_seurat,
    in  path image,
    in  string runtime,
    out map probabilities,
    out map annotated_seurat,
    out map batch_report,
    out map batch_csv,
    src py   "stages/assign_celltypes",
) split using (
    in  map annotated_seurat,
) using (
    threads = 30,
)
'''
def split(args):
    chunks = []
    for batch_id, path in args.annotated_seurat.items():
        chunk_def = {}
        chunk_def['annotated_seurat'] = path
        chunk_def['batch'] = batch_id
        chunk_def['__threads'] = 30
        chunk_def['__memgb'] = 16
        chunks.append(chunk_def)
    return {'chunks': chunks}

def format_cn_data(cn_file, filename):
    output = open(filename,"w")
    rows = open(cn_file,"r").read().splitlines()
    header = rows.pop(0).split("\t")
    discard = rows.pop(0)
    clones = header[4:]
    formatted_header = "chr,start,end,copy_number,clone"
    output.write(formatted_header+"\n")
    for row in rows:
        row = row.split("\t")
        prefix = row[:3]
        if prefix[1] == "-1": continue
        if prefix[0] == "X" or prefix[0] == "Y": continue
        prefix[0] = "chr" + prefix[0]
        for clone, cn in zip(clones, row[4:]):
            new_line = prefix + [cn, clone]
            output.write(",".join(new_line)+"\n")
    output.close()

def main(args, outs):
    # clone_annotated_seurat = "{}_clone_annotated.rds".format(args.batch)
    # clone_umap = "{}_clone_umap.svg".format(args.batch)
    # clonealign_fit = "{}_clone_umap.rds".format(args.batch)
    # clone_cn = "/work/shah/funnellt/projects/sc-mutsig/analysis/corrupt_tree/{0}/clones/{0}_padded_chrom_cn_clones.tsv".format(args.batch)
    # formatted_clone_cn = "{}_cn.csv".format(args.batch)
    # formatted_clone_cn = martian.make_path(formatted_clone_cn)
    # format_cn_data(clone_cn, formatted_clone_cn)
    # args.clone_cn = formatted_clone_cn
    # outs.clone_umap = martian.make_path(clone_umap)
    # outs.clone_annotated_seurat = martian.make_path(clone_annotated_seurat)
    # outs.clonealign_fit = martian.make_path(clonealign_fit)
    # scripts = scriptmanager.ScriptManager()
    # script = scripts.clonealign()
    # con = container.Container()
    # con.set_runtime(args.runtime)
    # con.set_image(args.image)
    # if os.path.exists("/work/shah/ceglian/pipeline_runs/test_clonealign/{}_clone_annotated.rds".format(args.batch)):
    #     outs.clone_annotated_seurat = "/work/shah/ceglian/pipeline_runs/test_clonealign/{}_clone_annotated.rds".format(args.batch)
    #     con.run(script, args, outs)
    pass
        

def join(args, outs, chunk_defs, chunk_outs):
    outs.clone_annotated_seurat = dict()
    outs.clonealign_fit = dict()
    outs.clone_umap = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        if out.clone_umap and os.path.exists(out.clone_umap):
            outs.clone_annotated_seurat[arg.batch] = out.clone_annotated_seurat
            outs.clone_umap[arg.batch] = out.clone_umap
            outs.clonealign_fit[arg.batch] = out.clonealign_fit
