import container
import scriptmanager
import martian

__MRO__ = '''
stage MERGE_SAMPLES(
    in  map qcd_scored_seurat,
    in  map batch,
    in  path image,
    in  string runtime,
    out map merged_seurat,
    out map merged_tsv,
    src py   "stages/merge_samples",
) split using (
    in  map batch,
) using (
    threads = 12,
)
'''
def split(args):
    chunks = []
    for batch_id, samples in args.batch.items():
        chunk_def = {}
        chunk_def['batched_samples'] = samples
        chunk_def['batch'] = batch_id
        chunk_def['__threads'] = 12
        chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    merged_batch = "{}_merged.rds".format(args.batch)
    merged_tsv = "{}_batch_samples.tsv".format(args.batch)
    outs.merged_seurat = martian.make_path(merged_batch)
    outs.merged_tsv = martian.make_path(merged_tsv)
    output = open(outs.merged_tsv,"w")
    output.write("sample_id\tpath\n")
    for sample, path in args.qcd_scored_seurat.items():
        if sample in args.batched_samples:
            output.write("{}\t{}\n".format(sample,path))
    output.close()
    scripts = scriptmanager.ScriptManager()
    script = scripts.mergesamples()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)

def join(args, outs, chunk_defs, chunk_outs):
    outs.merged_seurat = dict()
    outs.merged_tsv = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        outs.merged_seurat[arg.batch] = out.merged_seurat
        outs.merged_tsv[arg.batch] = out.merged_tsv
