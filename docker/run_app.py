import isabl_cli
import sys
import os
import subprocess
import argparse
import collections
import pprint

import metadata

pp = pprint.PrettyPrinter(indent=2)

CELLRANGER_VERSION = "3.1.0"
SCRNA_VERSION      = "2.0.0"

# todo add unit test
os.environ["ISABL_API_URL"]='https://isabl.shahlab.ca/api/v1/'
#os.environ["ISABL_API_URL"]='http://localhost:8000/api/v1/'
os.environ["ISABL_CLIENT_ID"]="1"

"""
#export ISABL_API_URL='https://isabl.shahlab.ca/api/v1/'
export ISABL_API_URL='http://localhost:8000/api/v1/'
export ISABL_CLIENT_ID=1
"""


def run(patient_id, run_cmd, app, app_status):
    valid_samples = metadata.patient_samples(patient_id)
    run_cmds = []
    for sample in valid_samples:
        unique_id = sample["unique_id"]
        if not app_status[unique_id]:
            experiment = isabl_cli.get_experiments(sample__identifier=unique_id)
            if len(experiment) > 0:
                systemid = experiment[0].system_id
                print("Found Experiment {} for {}".format(unique_id,systemid))
                cmd = run_cmd.format(system=systemid,app=app)
                run_cmds.append(cmd.split())
            else:
                print("No experiment found for {}".format(unique_id))
    if len(run_cmds) == 0:
        print("ALL SAMPLES FOR PATIENT {} HAVE BEEN RUN FOR {}.".format(patient_id,app))
    print("\n{} valid experiments found.".format(len(run_cmds)))
    print("\n\nStarting {} processes for isabl...".format(len(run_cmds)))
    for cmd in run_cmds:
        print("\n\t" + " ".join(cmd))
        subprocess.Popen(cmd)
    print("Completed.")

def check(patient_id, app, dependencies, force):
    '''
    Check and report status on console for all analyses related to specified patient, app and dependencies of the app.
    :param patient_id:
    :param app:
    :param dependencies:
    :return: a dictionary containing a list of app statuses, one items per sample {<sample_unique_id>: <message>, ...}
    '''
    app_status = dict()
    # get meta data for samples belonging to specified patient
    valid_samples = metadata.patient_samples(patient_id)
    # get all isabl analyses
    analyses = isabl_cli.get_analyses()
    for sample in valid_samples:
        applications = collections.defaultdict(int)
        dependent_status = collections.defaultdict(int)
        unique_id = sample["unique_id"]
        # get experiments for the specified sample
        experiments = isabl_cli.get_experiments(sample__identifier=unique_id)

        for analysis in analyses:
            if len(analysis.__dict__["targets"]) > 0:   # individual_centric and project centric analyses don't have targets (but, yes they should!)
                system_id = analysis.__dict__["targets"][0]["system_id"]
            else:
                continue
            for experiment in experiments:
                if system_id == experiment.system_id:
                    targetapp = analysis["application"].__dict__["name"]
                    targetversion = analysis["application"].__dict__["version"]
                    if app.upper() == targetapp and targetversion == CELLRANGER_VERSION:
                        status = analysis["status"]
                        applications[status] += 1
                    if app.upper() in dependencies and dependencies[app.upper()] == targetapp:
                        status = analysis["status"]
                        dependent_status[status] += 1

        if applications['SUCCEEDED'] == 0:
            print(unique_id, "no successes", applications)
            app_status[unique_id] = "**NOT COMPLETE**"
            if len(dependencies.keys()) > 0 and dependent_status["SUCCEEDED"] == 0:
                app_status[unique_id] = "**DEPENDENT APP NOT COMPLETE: {}**".format(dependencies[app.upper()])
        else:
            app_status[unique_id] = "SUCCEEDED"
        if applications['STARTED'] > 0:
            app_status[unique_id] = "*RUNNING*"
        if force:
            app_status[unique_id] = "**NOT COMPLETE**"
    return app_status

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--patient_id", type=str, help="Patient ID as found in google sheet (SPECTRUM-OV-XXX).")
    parser.add_argument("--app", type=str, help="Isabl app to run (cellranger or scrna).")
    parser.add_argument("--dependency", type=str, default=None, help="Isabl app dependencies for running app.")
    parser.add_argument("--status", action='store_true', help="Check and report status by sample only. Don't run apps.")
    parser.add_argument("--force", action='store_true', help="Force reanalysis")
    args = parser.parse_args()
    app = args.app.lower()
    dependencies = dict()

    if args.dependency:
        dependencies[app.upper()] = args.dependency.upper()

    assert app in ("scrna","cellranger"), "App must be 'scrna' or 'cellranger'."

    app_status = check(args.patient_id, app, dependencies, args.force)
    print("\n{app} STATUS BY SAMPLE:".format(app=app.upper()))
    for sample, msg in app_status.items():
        print(sample, app.upper(), msg)
        if "**NOT COMPLETE**" in msg:
            app_status[sample] = False
        else:
            app_status[sample] = True
    print("\n**********************")
    if not args.status:
        if app == "cellranger":
            app += "-{}".format(CELLRANGER_VERSION)
        elif app == 'scrna':
            app += "-{}".format(SCRNA_VERSION)
        print("\nRunning app {}".format(app))
        run_cmd = "isabl apps-grch38 {app} -fi system_id {system} --force"
        run(args.patient_id, run_cmd, app, app_status)
    else:
        print("Showing app status.")
