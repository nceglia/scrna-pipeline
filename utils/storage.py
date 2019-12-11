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

import os, uuid, sys
import subprocess
import json
import tarfile
import shutil
import gzip
from multiprocessing import Pool

def retrieve_cellranger_matrix(patient_id, standard_id):
    path_template = "{0}/outs/"
    valid_samples = metadata.patient_samples(patient_id)
    analyses = isabl_cli.get_analyses()
    analysis_stack = collections.defaultdict(list)
    for sample in valid_samples:
        if sample["nick_unique_id"] != standard_id: continue
        unique_id = sample["unique_id"]
        experiments = isabl_cli.get_experiments(sample__identifier=unique_id)
        for analysis in analyses:
            system_id = analysis.__dict__["targets"][0]["system_id"]
            for experiment in experiments:
                if system_id == experiment.system_id:
                    targetapp = analysis["application"].__dict__["name"]
                    if "CELLRANGER" == targetapp:
                        status = analysis["status"]
                        if status == "SUCCEEDED":
                            analysis_stack[unique_id].append((sample["nick_unique_id"], analysis["storage_url"]))
    return


class LocalTenxDataStorage(object):

    def __init__(self, sampleid):
        self.sampleid = sampleid.replace(".tar.gz","")
        if self.species == "mouse":
            self.sampleid += "_mouse"
        self.tenx_path = None
        self.cache = ".cache"
        try:
            os.makedirs(self.cache)
        except Exception as e:
            pass

    def download(self):



class TenxDataStorage(object):

    def __init__(self, sampleid, version="v3", software="cellranger", species="human"):
        self.species = species
        self.version = version
        self.sampleid = sampleid.replace(".tar.gz","")
        if self.species == "mouse":
            self.sampleid += "_mouse"
        self.storage_account = "scrnadata"
        if species != "human":
            self.container = "cellranger{}{}".format(species,version)
            self.bams = "{}bams".format(species)
            self.bus_container = "kallisto{}".format(species)
        else:
            self.container = "cellranger{}".format(version)
            self.bams = "bams"
            self.bus_container = "kallisto"
        self.block_blob_service = BlockBlobService(account_name='scrnadata', sas_token=aztok)
        self.tenx_path = None
        self.cache = ".cache"
        try:
            os.makedirs(self.cache)
        except Exception as e:
            pass

    def upload_cellranger(self, tenx):
        bam_tarball = tenx.bam_tarball()
        bam_tarball_name = os.path.split(bam_tarball)[-1]
        if self.sampleid not in bam_tarball:
            bam_tarball = os.path.join(self.sampleid, bam_tarball)
        print("Uploading {} {}".format(bam_tarball, bam_tarball_name))
        outs_tarball = tenx.outs_tarball()
        outs_tarball_name = os.path.split(outs_tarball)[-1]
        self.upload(self.bams, bam_tarball_name.replace("_mouse",""), bam_tarball)
        print(self.container, outs_tarball_name, outs_tarball)
        self.upload(self.container, outs_tarball_name.replace("_mouse",""), outs_tarball)

    def upload_kallisto(self, tenx):
        bus_tarball = tenx.bus_tarball()
        bus_tarball_name = os.path.split(bus_tarball)[-1]
        self.upload(self.bus_container, bus_tarball_name.replace("_mouse",""), bus_tarball)

    def unpack(self, path):
        tar = tarfile.open(path)
        tar.extractall(".cache/{}".format(self.sampleid))
        tar.close()
        if os.path.exists(".cache/outs"):
            os.rename(".cache/outs",os.path.join(self.cache,self.sampleid))

    def download(self):
        local = ".cache/{}.tar.gz".format(self.sampleid)
        local_unzipped = ".cache/{}".format(self.sampleid)
        gzipped = "{}.tar.gz".format(self.sampleid)
        print(gzipped, local, self.container)
        if not os.path.exists(local_unzipped):
            self.block_blob_service.get_blob_to_path(self.container, gzipped.replace("_mouse",""), local)
            self.unpack(local)
        self.tenx_path = os.path.join(self.cache,self.sampleid)
        return self.tenx_path

    def upload(self,container, blob, local):
        print ("Uploading {} to {} in {}".format(container,local,container))
        self.block_blob_service.create_blob_from_path(container, blob, local)
