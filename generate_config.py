import argparse
import yaml
import sys
import os

from utils.config import write_config

parser = argparse.ArgumentParser("SCRNA Pipeline")


parser.add_argument("--build", type=str, help="Reference build.", default="GRCh38")
parser.add_argument("--jobpath", type=str, help="Path to job folder.", default = "/results")
parser.add_argument("--aztok", type=str, help="Path to token folder", default = "/codebase/.sas_token")

subparsers = parser.add_subparsers(dest="subparser")

cellranger_args = subparsers.add_parser('cellranger')
cellranger_args.add_argument('--sampleid', type=str, help='Sample ID linked to fastqs in scrnadata.')
cellranger_args.add_argument("--datapath", type=str, help="Path to directory holding fastqs.", default ="/data")
cellranger_args.add_argument("--referencepath", type=str, help="Path to reference directory with reference in dir labelled by build.", default="/reference")
cellranger_args.add_argument("--cellranger", type=str, help="Path to cellranger bin.", default="/codebase/cellranger-3.0.2/cellranger-cs/3.0.2/bin/")

kallisto_args = subparsers.add_parser('kallisto')
kallisto_args.add_argument('--sampleid', type=str, help='Sample ID linked to fastqs in scrnadata.')
kallisto_args.add_argument("--datapath", type=str, help="Path to directory holding fastqs.", default ="/data")
kallisto_args.add_argument("--referencepath", type=str, help="Path to reference directory with reference in dir labelled by build.", default="/reference")
kallisto_args.add_argument("--kallisto", type=str, help="Path to kallisto bin.", default="/codebase/kallisto/")
kallisto_args.add_argument("--bustools", type=str, help="Path to bustools.", default="/codebase/bustools/")

sample_args = subparsers.add_parser('analysis')
sample_args.add_argument('--sampleid', type=str, help='Sample ID linked to fastqs in scrnadata.')
sample_args.add_argument("--markers", type=str, help="Path to marker matrix.")

clonealign_args = subparsers.add_parser('clonealign')
clonealign_args.add_argument('--samples', type=str, help='Sample ID linked to fastqs in scrnadata.')
clonealign_args.add_argument('--copy_cell_clones', type=str, help='Copy cell clones tsv.')
clonealign_args.add_argument('--filtered_cell_cn', type=str, help='Filter cell clones tsv.')
clonealign_args.add_argument("--markers", type=str, help="Path to marker matrix.")

integration_args = subparsers.add_parser("integrate")
integration_args.add_argument('--prefix', type=str, help='Analysis prefix')
integration_args.add_argument("--samples", type=str, help="txt file with valid samples on each line.")
integration_args.add_argument("--markers", type=str, help="Path to marker matrix.")
integration_args.add_argument("--datapath", type=str, help="Path to directory holding fastqs.", default ="/data")

integration_args = subparsers.add_parser("pipeline")
integration_args.add_argument('--sampleid', type=str, help='Sample ID linked to fastqs in scrnadata.')
integration_args.add_argument("--samples", type=str, help="txt file with valid samples on each line.", default="samples.txt")
integration_args.add_argument("--markers", type=str, help="Path to marker matrix.")
integration_args.add_argument("--datapath", type=str, help="Path to directory holding fastqs.", default ="/data")
integration_args.add_argument("--referencepath", type=str, help="Path to reference directory with reference in dir labelled by build.", default="/reference")
integration_args.add_argument("--cellranger", type=str, help="Path to cellranger bin.", default="/codebase/cellranger-3.0.2/cellranger-cs/3.0.2/bin/")

args = parser.parse_args()

if args.subparser != "clonealign" and args.subparser != "analysis" and not os.path.exists(args.samples):
    output = open(args.samples, "w")
    output.write(args.sampleid)
    output.close()

write_config(args)
print("Config written to settings.yaml")
