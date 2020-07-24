from jinja2 import Template
import subprocess
import pandas
import sys
import os

samples = eval(sys.argv[1])
project = sys.argv[2]
batches   = eval(sys.argv[3])
project_figure = sys.argv[4]
celltype_umaps = eval(sys.argv[5])
interactions = eval(sys.argv[6])
networks = eval(sys.argv[7])
markers = eval(sys.argv[8])
batch_markers = eval(sys.argv[9])
subtype_scores = eval(sys.argv[10])
enriched_pathways = eval(sys.argv[11])
ct_markers = sys.argv[12]
network_svg = eval(sys.argv[13])
frac_svg = eval(sys.argv[14])
copy_numbers = eval(sys.argv[15])
velocity = eval(sys.argv[16])
clone = eval(sys.argv[17])
report = sys.argv[18]

ctdf = pandas.read_csv(ct_markers)

sample_results = []
for sample in samples:
    svg_content = open(samples[sample],"r").read()
    if sample in frac_svg:
        frac = open(frac_svg[sample],"r").read()
    else:
        frac = ""
    sample_results.append({"name":sample,"svg":svg_content,"frac":frac})
batch_results = []
for batch in batches:
    svg_content = open(batches[batch],"r").read()
    bm = batch_markers[batch]
    df = pandas.read_csv(bm)
    if batch in copy_numbers:
        subprocess.check_output(["pdf2svg",copy_numbers[batch], copy_numbers[batch].replace("pdf","svg")])
        copy_number_svg = open(copy_numbers[batch].replace("pdf","svg"),"r").read()
    else:
        copy_number_svg = ""
    if batch in clone:
        clone_svg = open(clone[batch],"r").read()
    else:
        clone_svg = ""
    batch_results.append({"name":batch, "svg":svg_content, "copy_number": copy_number_svg, "table": df.to_html(border=0, table_id="batch_markers_{}".format(batch), classes="batch_markers"), "id": "batch_markers_{}".format(batch), "pdf": copy_numbers[batch], "clone": clone_svg})

celltype_results = []
for celltype in celltype_umaps:
    svg_content = open(celltype_umaps[celltype],"r").read()
    velocity_svg = open(velocity[celltype]["svg"],"r").read()
    bm = markers[celltype]
    df = pandas.read_csv(bm)
    clusters = []
    if celltype in enriched_pathways:
        for cluster in enriched_pathways[celltype]:
            svg_content_pathway = open(cluster["svg"],"r").read()
            clusters.append({"cluster":cluster["cluster"],"svg":svg_content_pathway})
    celltype = celltype.replace(".","-").capitalize()
    celltype_results.append({"name":celltype, "svg":svg_content, "table": df.to_html(border=0, classes="batch_markers", table_id="celltype_markers_{}".format(celltype)), "id":"celltype_markers_{}".format(celltype), "clusters": clusters, "velocity": velocity_svg})
    

interaction_results = []
for pathway, batches in interactions.items():
    results = []
    for batch, result in batches.items():
        svg_content_pathway = open(result["svg"],"r").read()
        svg_content_net = open(network_svg[pathway][batch]["svg"],"r").read()
        results.append({"batch":batch, "svg":svg_content_pathway, "svg_net": svg_content_net})
    interaction_results.append({"pathway":pathway,"results":results})

subtype_results = []
for subtype, result in subtype_scores.items():
    svg_content = open(result["svg"],"r").read()
    celltype, subtype = subtype.split("_")
    subtype = "{}-{}".format(subtype, celltype.replace(".","-"))
    subtype_results.append({"name":subtype,"svg":svg_content})

svg_project_figure = open(project_figure,"r").read()
project_markers = ctdf.to_html(border=0, table_id="batch_markers_project", classes="batch_markers")

cwd = os.path.dirname(os.path.abspath(__file__))
template = Template(open(os.path.join(cwd,"../html/report.html"),"r").read())

html = template.render(samples=sample_results,
                       project=project,
                       batches=batch_results,
                       project_figure=svg_project_figure,
                       celltypes=celltype_results,
                       pathways=interaction_results,
                       project_markers=project_markers,
                       subtypes=subtype_results)
output = open(report,"w")
output.write(html)
output.close()
