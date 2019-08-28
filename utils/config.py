import os
import yaml

from utils.cloud import ReferenceDataStorage

yaml_file = os.path.join(os.getcwd(), "settings.yaml")
def yaml_configuration():
    if os.path.exists(yaml_file):
        with open(yaml_file, "r") as f:
            doc = yaml.load(f)
            return doc

basic_yaml = """
run_command: {run_command}
prefix: "{prefix}"
samples: "{samples}"
build: "{build}"
jobpath: "{jobpath}"
datapath: "{datapath}"
referencepath: "{referencepath}"
rho_matrix: {markers}
cellranger: "{cellranger}"
copy_number_data: None
scviz_embedding: None
run_scvis: False
run_cellassign: False
run_clonealign: False
plot_scvis: False
clustering: False
report: True
perplexity: 5
resolution: 0.2
stds: 6
components: 50
chemistry: auto
low_counts_genes_threshold: 4
qc_type: "standard"
scviz_embedding: None
copy_number_data: None
mito: 20
mouse_reference: /reference/mm10
kallisto: {kallisto}
bustools: {bustools}
copy_cell_clones: {copy_cell_clones}
filtered_cell_cn: {filtered_cell_cn}
"""


def write_config(args):
    output = open("settings.yaml","w")
    subcommand = args.subparser
    jobpath    = args.jobpath
    build = args.build
    prefix = "scrna"
    cellranger = "none"
    referencepath = "none"
    markers = "none"
    datapath = "none"
    samples = "none"
    sample_file = "none"
    datapath = "/data"
    kallisto = "/codebase/kallisto"
    bustools = "/codebase/bustools"
    filtered_cell_cn = "none"
    copy_cell_clones = "none"

    if subcommand == "cellranger":
        datapath = args.datapath
        referencepath = args.referencepath
        cellranger = args.cellranger
        prefix     = args.sampleid
        datapath = args.datapath
    elif subcommand == "sample":
        markers = args.markers
        prefix  = args.sampleid
    elif subcommand == "integrate":
        sample_file = args.samples
        prefix = args.prefix
        markers = args.markers
    elif subcommand == "pipeline":
        sample_file = args.samples
        prefix = args.sampleid
        markers = args.markers
        datapath = args.datapath
        referencepath = args.referencepath
        cellranger = args.cellranger
    elif subcommand == "kallisto":
        sample_file = args.samples
        prefix = args.sampleid
        markers = args.markers
        datapath = args.datapath
        referencepath = args.referencepath
        kallisto = args.kallisto
        bustools = args.bustools
    elif subcommand == "clonealign":
        prefix = args.sampleid
        markers = args.markers
        copy_cell_clones = args.copy_cell_clones
        filtered_cell_cn = args.filtered_cell_cn

    output.write(basic_yaml.format(run_command=subcommand,
                                   prefix=prefix,
                                   build=build,
                                   jobpath=jobpath,
                                   datapath=datapath,
                                   referencepath=referencepath,
                                   cellranger=cellranger,
                                   kallisto=kallisto,
                                   markers=markers,
                                   samples=sample_file,
                                   bustools=bustools,
                                   copy_cell_clones=copy_cell_clones,
                                   filtered_cell_cn=filtered_cell_cn))

class Configuration(object):
    def __init__(self):
        self.build = "GRCh38"
        overrides = yaml_configuration()
        if overrides != None:
            for attr, value in overrides.items():
                setattr(self, attr, value)
        if self.run_command == "cellranger" or self.run_command == "pipeline":
            refobj = ReferenceDataStorage(self.build, self.referencepath)
            if not hasattr(self, "reference"):
                self.reference = refobj.download()
            self.genes_gtf = os.path.join(self.reference, "genes/genes.gtf")
