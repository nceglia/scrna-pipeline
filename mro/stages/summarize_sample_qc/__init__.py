import container
import scriptmanager
import martian
import itertools

__MRO__ = '''
stage SUMMARIZE_SAMPLE_QC(
    in  map qcd_scored_seurat,
    in  map csv,
    in  map sce,
    in  path markers,
    in  path image,
    in  string runtime,
    out png qc_report,
    src py "stages/summarize_sample_qc",
) split using (
    in  map qcd_scored_seurat,
) using (
    threads = 8,
)
'''

def split(args):
    chunks = []
    for sample, path in args.qcd_scored_seurat.items():
        chunk_def = {}
        chunk_def["seurat"]   = path
        chunk_def["sce"]      = args.sce[sample]
        chunk_def["csv"]      = args.csv[sample]
        chunk_def["markers"]  = args.markers
        chunk_def["sample"]   = sample
        chunk_def['__threads'] = 8
        chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    png = "{}_qc_report.svg".format(args.sample)
    outs.qc_report = martian.make_path(png)
    scripts = scriptmanager.ScriptManager()
    script = scripts.summarizeqc()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)

def join(args, outs, chunk_defs, chunk_outs):
    # cwd = os.path.dirname(os.path.abspath(__file__))
    # template = Template(open(os.path.join(cwd,"../../../lib/html/report.html"),"r").read())
    # sample_results = []
    outs.qc_report = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        outs.qc_report[arg.sample] = out.qc_report
    # html = template.render(samples=sample_results)
    # output = open(outs.qc_report,"w")
    # output.write(html)
    # output.close()
