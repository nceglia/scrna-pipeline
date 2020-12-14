import container
import scriptmanager
import martian
import glob
import os
import shutil

__MRO__ = '''
stage BATCH_CORRECTION(
    in  rds batch_merged_seurat,
    in  path image,
    in  string runtime,
    out rds integrated_seurat,
    out svg project_figure,
    src py   "stages/batch_correction",
) using (
    threads = 16,
)
'''

def main(args, outs):
    scripts = scriptmanager.ScriptManager()
    script = scripts.batchcorrection()
    con = container.Container()
    cell_annotations = "cells.tsv"
    outs.cell_annotations = martian.make_path(cell_annotations)
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)
    base = outs.cell_annotations.split("RNASCP")[0]
    copied = os.path.join(base, cell_annotations)
    shutil.copyfile(outs.cell_annotations, copied)