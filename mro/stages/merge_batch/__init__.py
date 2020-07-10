import container
import scriptmanager
import martian

__MRO__ = '''
stage MERGE_BATCHES(
    in  map annotated_seurat,
    in  path image,
    in  string runtime,
    out rds batch_merged_seurat,
    out path batch_tsv,
    out path celltype_csv,
    out string[] celltypes,
    src py   "stages/merge_batch",
) using (
    threads = 16,
)
'''

def main(args, outs):
    output = open(outs.batch_tsv,"w")
    output.write("batch_id\tpath\n")
    for batch_id, path in args.annotated_seurat.items():
        output.write("{}\t{}\n".format(batch_id,path))
    output.close()
    outs.batch_merged_seurat = "/juno/work/shah/ceglian/rnascp/signatures_seurats/batch_merged_seurat.rds"
    scripts = scriptmanager.ScriptManager()
    script = scripts.mergebatches()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)
    celltypes = open(outs.celltype_csv,"r").read().splitlines()
    celltypes.pop(0)
    valid_celltypes = []
    for celltype in celltypes:
        ctype = celltype.split(",")[1].replace('"',"")
        count = int(celltype.split(",")[2])
        if count > 50:
            valid_celltypes.append(ctype)
    outs.celltypes = valid_celltypes
