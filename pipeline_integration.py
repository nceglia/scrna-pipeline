import argparse
import glob
import os
import yaml
import sys

import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

from utils.config import Configuration
from utils.cloud import SampleCollection

from tasks.run_cellranger import RunCellranger
from tasks.run_qc import RunQC
from tasks.run_report import RunReport
from tasks.run_cellassign import RunCellAssign
from tasks.run_clustering import RunClustering
from tasks.run_collection import RunCollection

config = Configuration()

def create_workflow():
    print("Creating integration analysis workflow.")
    workflow = pypeliner.workflow.Workflow()

    prefix = config.prefix
    print("******* {} ********".format(prefix))

    workflow = RunCollection(workflow)
    # workflow = RunCellAssign(prefix, workflow)
    # workflow = RunReport(prefix, workflow)

    return workflow

if __name__ == '__main__':

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    parsed_args = argparser.parse_args()

    args = vars(parsed_args)
    workflow = create_workflow()
    pyp = pypeliner.app.Pypeline(config=args)
    pyp.run(workflow)
