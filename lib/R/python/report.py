from jinja2 import Template
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
report = sys.argv[13]

ctdf = pandas.read_csv(ct_markers)

sample_results = []
for sample in samples:
    svg_content = open(samples[sample],"r").read()
    sample_results.append({"name":sample,"svg":svg_content})

batch_results = []
for batch in batches:
    svg_content = open(batches[batch],"r").read()
    bm = batch_markers[batch]
    df = pandas.read_csv(bm)
    batch_results.append({"name":batch, "svg":svg_content, "table": df.to_html(border=0, table_id="batch_markers_{}".format(batch), classes="batch_markers"), "id": "batch_markers_{}".format(batch)})

celltype_results = []
for celltype in celltype_umaps:
    svg_content = open(celltype_umaps[celltype],"r").read()
    bm = markers[celltype]
    df = pandas.read_csv(bm)
    if celltype in networks:
        svg_content_network = open(networks[celltype],"r").read()
    else:
        svg_content_network = ""
    if celltype in enriched_pathways:
        svg_content_pathway = open(enriched_pathways[celltype],"r").read()
    else:
        svg_content_pathway = ""
    celltype = celltype.replace(".","-").capitalize()
    celltype_results.append({"name":celltype, "svg":svg_content, "table": df.to_html(border=0, classes="batch_markers", table_id="celltype_markers_{}".format(celltype)), "id":"celltype_markers_{}".format(celltype), "enrichment_svg": svg_content_network, "pathway_svg": svg_content_pathway})

interaction_results = []
for pathway, results in interactions.items():
    svg_contents = []
    for result in results:
        svg_contents.append(open(result["svg"],"r").read())
    interaction_results.append({"pathway":pathway,"svgs":svg_contents})

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
