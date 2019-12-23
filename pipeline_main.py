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

from workflows.run_main import RunMain

config = Configuration()

def create_workflow():
    print("Creating integration analysis workflow.")
    workflow = pypeliner.workflow.Workflow()

    prefix = config.prefix
    filtered_matrices = config.matrix
    print("******* {} ********".format(prefix))
    print("*** {}".format(filtered_matrices))
    workflow = RunMain(workflow)

    return workflow

if __name__ == '__main__':

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    parsed_args = argparser.parse_args()

    args = vars(parsed_args)
    workflow = create_workflow()
    pyp = pypeliner.app.Pypeline(config=args)
    pyp.run(workflow)
