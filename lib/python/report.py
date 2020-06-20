from jinja2 import Template
import sys
import os

samples = eval(sys.argv[1])
project = sys.argv[2]
batches   = eval(sys.argv[3])
sample_composition = sys.argv[4]
batch_composition = sys.argv[5]
umap_integrated = sys.argv[6]
umap_merged = sys.argv[7]
celltype_umaps = eval(sys.argv[8])
interactions = eval(sys.argv[9])
networks = eval(sys.argv[10])
markers = eval(sys.argv[11])
report = sys.argv[12]

sample_results = []
for sample in samples:
    svg_content = open(samples[sample],"r").read()
    sample_results.append({"name":sample,"svg":svg_content})

batch_results = []
for batch in batches:
    svg_content = open(batches[batch],"r").read()
    batch_results.append({"name":batch, "svg":svg_content})

celltype_results = []
for celltype in celltype_umaps:
    svg_content = open(celltype_umaps[celltype],"r").read()
    celltype_results.append({"name":celltype, "svg":svg_content})

network_results = []
for celltype in networks:
    for cluster in networks[celltype]:
        svg_content = open(networks[celltype][cluster],"r").read()
        svg_content = '<svg' + svg_content.split('<svg')[1]
        network_results.append({"name":celltype,"cluster":cluster,"svg":svg_content})

interaction_results = []
for sender in interactions:
    for receiver in interactions[sender]:
        for pathway in interactions[sender][receiver]:
            svg_content = open(interactions[sender][receiver][pathway],"r").read()
            interaction_results.append({"sender":sender,"receiver":receiver,"pathway":pathway,"svg":svg_content})

svg_sample_composition = open(sample_composition,"r").read()
svg_batch_composition = open(batch_composition,"r").read()
svg_umap_integrated = open(umap_integrated,"r").read()
svg_umap_merged = open(umap_merged,"r").read()

cwd = os.path.dirname(os.path.abspath(__file__))
template = Template(open(os.path.join(cwd,"../html/report.html"),"r").read())

html = template.render(samples=sample_results,
                       project=project,
                       sample_composition=svg_sample_composition,
                       batches=batch_results,
                       umap_integrated=svg_umap_integrated,
                       umap_merged=svg_umap_merged,
                       batch_composition=svg_batch_composition,
                       celltypes=celltype_results,
                       interactions=interaction_results,
                       networks=network_results)
output = open(report,"w")
output.write(html)
output.close()
