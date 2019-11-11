import subprocess
import glob
import collections
import argparse
import metadata
import isabl_cli
import os
import gzip
import shutil
import pickle
import metadata


class TenxDataStorage(object):

    def __init__(self, sampleid):
        self.sampleid = sampleid.replace(".tar.gz","")
        self.tenx_path = None
        # self.cache = ".cache"
        # try:
        #     os.makedirs(self.cache)
        # except Exception as e:
        #     pass

    def retrieve_latest(self):
        path_template = "{0}/outs/"
        sample = metadata.single_sample(self.sampleid)
        analyses = isabl_cli.get_analyses()
        analysis_stack = collections.defaultdict(list)

        unique_id = sample["unique_id"]
        experiments = isabl_cli.get_experiments(sample__identifier=unique_id)
        for analysis in analyses:
            system_id = analysis.__dict__["targets"][0]["system_id"]
            for experiment in experiments:
                if system_id != experiment.system_id:
                    continue
                if "CELLRANGER" != analysis["application"].__dict__["name"]:
                    continue
                print("Cellranger")
                print(analysis["status"])
                if analysis["status"] != "SUCCEEDED":
                    print("sucess")
                    analysis_stack[unique_id].append(analysis["storage_url"])
                    print("hit")

        paths = dict()
        for sample, urls in analysis_stack.items():
            latest = urls.pop()
            cellranger_outs = os.path.join(latest, path_template.format(sample))
            if os.path.exists(cellranger_outs):
                print(cellranger_outs)
                break

        # return paths
        #
        # self.tenx_path = os.path.join(self.cache,self.sampleid)
        # return self.tenx_path


if __name__ == '__main__':
    storage = TenxDataStorage("SPECTRUM-OV-052_S1_CD45N_INFRACOLIC_OMENTUM")
    storage.retrieve_latest()
