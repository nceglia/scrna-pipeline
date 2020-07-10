import container
import scriptmanager
import martian

__MRO__ = '''
stage GENERATE_REPORT(
    in  string project,
    in  map qc_report,
    in  map batch_report,
    in  svg sample_composition,
    in  svg batch_composition,
    in  svg umap_integrated,
    in  svg umap_merged,
    in  map celltype_umaps,
    in  map interactions,
    in  map networks,
    in  map markers,
    in  path image,
    in  string runtime,
    out html report,
    src py   "stages/generate_report",
) using (
    threads = 12,
)
'''

def main(args, outs):
    scripts = scriptmanager.ScriptManager()
    script = scripts.generatereport()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)
    raise ValueError(outs)