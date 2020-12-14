import container
import scriptmanager
import martian
import os

__MRO__ = '''
stage MERGE_BATCHES(
    in  map annotated_sce,
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
    for batch_id, path in args.annotated_sce.items():
        if os.path.exists(path):
            output.write("{}\t{}\n".format(batch_id,path))
    output.close()
    scripts = scriptmanager.ScriptManager()
    script = scripts.mergebatches()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)
