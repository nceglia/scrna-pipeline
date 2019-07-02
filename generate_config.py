import argparse
import yaml
import sys
import os

from utils.config import write_config

parser = argparse.ArgumentParser("SCRNA Pipeline")


parser.add_argument("--build", type=str, help="Reference build.", default="GRCh38")
parser.add_argument("--jobpath", type=str, help="Path to job folder.", default = "/results")

subparsers = parser.add_subparsers(dest="subparser")

cellranger_args = subparsers.add_parser('cellranger')
cellranger_args.add_argument('--sampleid', type=str, help='Sample ID linked to fastqs in scrnadata.')
cellranger_args.add_argument("--datapath", type=str, help="Path to directory holding fastqs.", default ="/data")
cellranger_args.add_argument("--referencepath", type=str, help="Path to reference directory with reference in dir labelled by build.", default="/reference")
cellranger_args.add_argument("--cellranger", type=str, help="Path to cellranger bin.", default="/codebase/cellranger-3.0.2/cellranger-cs/3.0.2/bin/")

sample_args = subparsers.add_parser('analysis')
sample_args.add_argument('--sampleid', type=str, help='Sample ID linked to fastqs in scrnadata.')
sample_args.add_argument("--markers", type=str, help="Path to marker matrix.")

integration_args = subparsers.add_parser("integrate")
integration_args.add_argument('--prefix', type=str, help='Analysis prefix')
integration_args.add_argument("--samples", type=str, help="txt file with valid samples on each line.")
integration_args.add_argument("--markers", type=str, help="Path to marker matrix.")

args = parser.parse_args()

write_config(args)
print("Config written to settings.yaml")
