from azure.common.client_factory import get_client_from_cli_profile
from azure.mgmt.compute import ComputeManagementClient
from azure.mgmt.network import NetworkManagementClient
import os, uuid, sys
from azure.storage.blob import BlockBlobService, PublicAccess
import subprocess
import json
import tarfile
import shutil
import gzip
from multiprocessing import Pool


def get_tenx_object(sampleid):
    tenx = TenxDataStorage(sampleid)
    path = tenx.download()
    return path

aztok = open(".sas_token","r").read().strip()

class SampleCollection(object):

    def __init__(self, sampleids):
        self.tenxs = sampleids
        pool = Pool(processes=len(self.tenxs))
        paths = pool.map(get_tenx_object, self.tenxs)
        self.samples = dict(zip(self.tenxs,paths))
        print("Finished cache")

    def sample_paths(self):
        return self.samples


class TenxDataStorage(object):

    def __init__(self, sampleid, version="v3"):
        self.version = version
        self.sampleid = sampleid.replace(".tar.gz","")
        self.storage_account = "scrnadata"
        self.container = "cellranger{}".format(version)
        self.rdatacontainer = "rawrdata{}".format(version)
        self.block_blob_service = BlockBlobService(account_name='scrnadata', sas_token=aztok)
        self.tenx_path = None
        self.cache = ".cache"
        try:
            os.makedirs(self.cache)
        except Exception as e:
            pass

    def rdata(self):
        local = ".cache/{}.rdata".format(self.sampleid)
        raw = "{}.rdata".format(self.sampleid)
        if not os.path.exists(local):
            self.block_blob_service.get_blob_to_path(self.rdatacontainer, raw, local)
        return local

    def upload_cellranger(self, tenx):
        bam_tarball = tenx.bam_tarball()
        bam_tarball_name = os.path.split(bam_tarball)[-1]
        outs_tarball = tenx.outs_tarball()
        outs_tarball_name = os.path.split(outs_tarball)[-1]
        self.upload("bams" ,bam_tarball_name, bam_tarball)
        self.upload("cellranger{}".format(self.version), outs_tarball_name, outs_tarball)


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
        if not os.path.exists(local_unzipped):
            self.block_blob_service.get_blob_to_path(self.container, gzipped, local)
            self.unpack(local)
        self.tenx_path = os.path.join(self.cache,self.sampleid)
        return self.tenx_path

    def upload(self,container, blob, local):
        print ("Uploading {} to {} in {}".format(container,local,container))
        self.block_blob_service.create_blob_from_path(container, blob, local)


class RdataStorage(object):

    def __init__(self, cache=".cache/sces"):
        if not os.path.exists(cache):
            os.makedirs(cache)
        self.cache = cache
        self.storage_account = "scrnadata"
        self.container = "rdatarawv3"
        self.block_blob_service = BlockBlobService(account_name=self.storage_account, sas_token=aztok)

    def get_sce(self,  sampleid):
        rdata = "{}.rdata".format(sampleid)
        local = os.path.join(self.cache,rdata)
        self.block_blob_service.get_blob_to_path(self.container, rdata, local)
        assert os.path.exists(local), "Download didn't work."
        return local

class ReportStorage(object):

    def __init__(self, results_path):
        self.results_path = results_path
        self.storage_account = "scrnadata"
        self.container = "reports"
        self.block_blob_service = BlockBlobService(account_name=self.storage_account, sas_token=aztok)

    def upload(self, path, sampleid):
        reporttarball = os.path.join(path, "{}.tar.gz".format(sampleid))
        with tarfile.open(reporttarball, "w:gz") as tar:
            tar.add(self.results_path, arcname=os.path.basename(reporttarball))
        blobtarball = os.path.split(reporttarball)[1]
        self.block_blob_service.create_blob_from_path(self.container, blobtarball, reporttarball)

class FastqDataStorage(object):

    def __init__(self, sampleid):
        self.sampleid = sampleid
        self.storage_account = "scrnadata"
        self.container = "rnaseq"
        self.block_blob_service = BlockBlobService(account_name=self.storage_account, sas_token=aztok)

    def set_data_path(self, path):
        self.datapath = path
        try:
            os.makedirs(path)
        except Exception as e:
            pass

    def download_fastqs(self):
        fastqs = self.block_blob_service.list_blobs(self.container)
        for fastq in fastqs:
            if len(fastq.name.split("/")) != 2: continue
            sample, fastqname = fastq.name.split("/")
            if sample == self.sampleid:
                local = os.path.join(self.datapath, fastqname)
                if not os.path.exists(local):
                    self.block_blob_service.get_blob_to_path(self.container, fastq.name, local)
        return self.datapath

class ReferenceDataStorage(object):

    def __init__(self,build, referencepath):
        self.sampleid = build
        self.storage_account = "scrnadata"
        self.container = "reference"
        self.block_blob_service = BlockBlobService(account_name='scrnadata', sas_token=aztok)
        self.reference = os.path.join(referencepath, build)
        self.build = build
        self.referencepath = referencepath

    def extract(self,path):
        tar = tarfile.open(path)
        tar.extractall(path=self.referencepath)
        tar.close()

    def download(self):
        local = os.path.join(self.referencepath, "{}.tar.gz".format(self.build))
        if not os.path.exists(self.reference):
            self.block_blob_service.get_blob_to_path(self.container, "{}.tar.gz".format(self.build), local)
            self.extract(local)
        return self.reference


class VirtualMachine(object):

    def __init__(self):
        compute_client = get_client_from_cli_profile(ComputeManagementClient)
        network_client = get_client_from_cli_profile(NetworkManagementClient)
        nic = network_client.network_interfaces.get(GROUP_NAME, NIC_NAME)
        print('\nCreating Linux Virtual Machine')
        vm_parameters = create_vm_parameters(nic.id, VM_REFERENCE['linux'])
        print(vm_parameters)
        async_vm_creation = compute_client.virtual_machines.create_or_update(GROUP_NAME, VM_NAME, vm_parameters)
        async_vm_creation.wait()

    def run_job(self):
        pass

    def check_status(self):
        pass

if __name__ == '__main__':
    test_fastqs = FastqDataStorage("test-fastqs")
    test_fastqs.download_fastqs()
