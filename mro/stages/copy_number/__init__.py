import container
import scriptmanager
import martian
import os
import shutil

__MRO__ = '''
stage COPY_NUMBER(
    in  map annotated_seurat,
    in  path cn_ref_genes,
    in  path cn_reference,
    in  path image,
    in  string runtime,
    out map cn_output,
    src py   "stages/copy_number",
) split using (
    in  map annotated_seurat,
) using (
    threads = 30,
)
'''

def split(args):
    chunks = []
    if not args.cn_ref_genes or not args.cn_ref_genes:
        return {"chunks": chunks}
    for batch, obj in args.annotated_seurat.items():
        chunk_def = {}
        chunk_def['batch'] = batch
        chunk_def['seurat'] = obj
        chunk_def['__threads'] = 30
        chunk_def['__memgb'] = 12
        chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    args.cn_reference = "/juno/work/shah/ceglian/rnascp/resources/spectrum_infercnv_reference_20200526.rdata"
    path = "{}_infercnv".format(args.batch)
    outs.cn_output = martian.make_path(path)
    scripts = scriptmanager.ScriptManager()
    script = scripts.infercnv()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    # con.run(script, args, outs)

def join(args, outs, chunk_defs, chunk_outs):
    outs.cn_output = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        pdf = os.path.join(out.cn_output,"infercnv.pdf")
        if os.path.exists(pdf):
            outs.cn_output[arg.batch] = pdf
