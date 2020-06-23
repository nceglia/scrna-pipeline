import container
import scriptmanager
import martian

__MRO__ = '''
stage ASSIGN_CELLTYPES(
    in  map merged_seurat,
    in  path marker_csv,
    in  path image,
    in  string runtime,
    out map probabilities,
    out map annotated_seurat,
    out map batch_report,
    out map batch_csv,
    src py   "stages/assign_celltypes",
) split using (
    in  map merged_seurat,
) using (
    threads = 30,
)
'''
def split(args):
    chunks = []
    for batch_id, path in args.merged_seurat.items():
        chunk_def = {}
        chunk_def['merged_seurat'] = path
        chunk_def['batch'] = batch_id
        chunk_def['__threads'] = 30
        chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    annotated_seurat = "{}_annotated.rds".format(args.batch)
    batch_report = "{}_report.svg".format(args.batch)
    probabilities = "{}_probabilities.rds".format(args.batch)
    batch_markers = "{}_markers.csv".format(outs.batch_csv)
    outs.batch_report = martian.make_path(batch_report)
    outs.batch_csv = martian.make_path(batch_markers)
    outs.annotated_seurat = martian.make_path(annotated_seurat)
    outs.probabilities = martian.make_path(probabilities)
    scripts = scriptmanager.ScriptManager()
    script = scripts.assigncelltypes()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)

def join(args, outs, chunk_defs, chunk_outs):
    outs.probabilities = dict()
    outs.annotated_seurat = dict()
    outs.batch_report = dict()
    outs.batch_csv = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        outs.annotated_seurat[arg.batch] = out.annotated_seurat
        outs.probabilities[arg.batch] = out.probabilities
        outs.batch_csv[arg.batch] = out.batch_csv
        outs.batch_report[arg.batch] = out.batch_report
