import container
import scriptmanager
import martian
import os

__MRO__ = '''
stage ENRICHMENT_NETWORK(
    in  map markers,
    in  path gmt,
    in  path image,
    in  string runtime,
    out map pathway_network,
    src py   "stages/enrichment_network",
) split using (
    in  map markers,
) using (
    threads = 4,
)
'''
#"","p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","gene"
def clusters(csv):
    lines = open(csv,"r").read().splitlines()
    header = lines.pop(0)
    _clusters = []
    for line in lines:
        line = line.split(",")
        cluster = line[-2].replace('"',"")
        _clusters.append(cluster)
    return set(_clusters)


def split(args):
    chunks = []
    cluster_labels = clusters(args.ct_markers)
    for cluster in cluster_labels:
        chunk_def = {}
        chunk_def['markers'] = args.ct_markers
        chunk_def['celltype'] = cluster
        chunk_def['__threads'] = 4
        chunks.append(chunk_def)
    return {'chunks': chunks}

def main(args, outs):
    png = "{}_network.svg".format(args.celltype)
    outs.pathway_network = martian.make_path(png)
    png = "{}_enriched.svg".format(args.celltype)
    outs.enriched_pathways = martian.make_path(png)
    scripts = scriptmanager.ScriptManager()
    script = scripts.enrichmentnetwork()
    con = container.Container()
    con.set_runtime(args.runtime)
    con.set_image(args.image)
    con.run(script, args, outs)

def join(args, outs, chunk_defs, chunk_outs):
    outs.pathway_network = dict()
    outs.enriched_pathways = dict()
    for arg, out in zip(chunk_defs, chunk_outs):
        if os.path.exists(out.pathway_network) and not os.stat(out.pathway_network).st_size == 0:
            outs.pathway_network[arg.celltype] = out.pathway_network
            outs.enriched_pathways[arg.celltype] = out.enriched_pathways
