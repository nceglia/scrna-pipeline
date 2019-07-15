from azure.common.client_factory import get_client_from_cli_profile
from azure.mgmt.compute import ComputeManagementClient
from azure.mgmt.network import NetworkManagementClient
from azure.common.credentials import ServicePrincipalCredentials
from azure.mgmt.resource import ResourceManagementClient
from azure.mgmt.compute.models import DiskCreateOption

from msrestazure.azure_exceptions import CloudError


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

    def __init__(self, sampleid, version="v3", species="human"):
        self.species = species
        self.version = version
        self.sampleid = sampleid.replace(".tar.gz","")
        if self.species == "mouse":
            self.sampleid += "_mouse"
        self.storage_account = "scrnadata"
        if species != "human":
            self.container = "cellranger{}{}".format(species,version)
            self.bams = "{}bams".format(species)
            self.rdatacontainer = "rdataraw{}{}".format(species,version)
        else:
            self.container = "cellranger{}".format(version)
            self.bams = "bams"
            self.rdatacontainer = "rdataraw{}".format(version)
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
        # self.upload(self.bams, bam_tarball_name, bam_tarball)
        # self.upload(self.container, outs_tarball_name, outs_tarball)


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
        uploaded = False
        for blobdef in self.block_blob_service.list_blobs(self.container):
            if blob in blodef.name:
                uploaded = True
        if not uploaded:
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
        self.LOCATION = 'canadacentral'
        self.GROUP_NAME = 'nick'
        self.VNET_NAME = 'nick-vnet'
        self.NIC_NAME = 'nickvm405'
        self.VM_NAME = 'scrna-pipeline3'
        self.USERNAME = 'scrna-user'
        self.PASSWORD = 'Password123!'
        self.IP_ADDRESS_NAME='nickvm-ip'
        self.VM_REFERENCE = {
            'linux': {
                'publisher': 'Canonical',
                'offer': 'UbuntuServer',
                'sku': '16.04.0-LTS',
                'version': 'latest'
            },
        }
        self.compute_client = get_client_from_cli_profile(ComputeManagementClient)
        self.network_client = get_client_from_cli_profile(NetworkManagementClient)
        self.nic = self.network_client.network_interfaces.get(self.GROUP_NAME, self.NIC_NAME)
        print('\nCreating Linux Virtual Machine')
        vm_parameters = self.create_vm_parameters()
        async_vm_creation = self.compute_client.virtual_machines.create_or_update(self.GROUP_NAME, self.VM_NAME, vm_parameters)
        print("\nSpinning up FTL.")
        async_vm_creation.wait()


    def create_vm_parameters(self):
        vm_reference = self.VM_REFERENCE["linux"]
        return {
            'location': self.LOCATION,
            'os_profile': {
                'computer_name': self.VM_NAME,
                'admin_username': self.USERNAME,
                'admin_password': self.PASSWORD
            },
            'hardware_profile': {
                'vm_size': 'Standard_DS1_v2'
            },
            'storage_profile': {
                'image_reference': {
                    'publisher': vm_reference['publisher'],
                    'offer': vm_reference['offer'],
                    'sku': vm_reference['sku'],
                    'version': vm_reference['version']
                },
            },
            'network_profile': {
                'network_interfaces': [{
                    'id': self.nic.id,
                }]
            },
        }

    def start(self):
        print("Warming up...")
        async_vm_start = self.compute_client.virtual_machines.start(self.GROUP_NAME, self.VM_NAME)
        async_vm_start.wait()

    def deallocate(self):
        print("Cooling down...")
        async_vm_deallocate = self.compute_client.virtual_machines.deallocate(self.GROUP_NAME, self.VM_NAME)
        async_vm_deallocate.wait()

    def stop(self):
        print("Stopping engine.")
        async_vm_stop = self.compute_client.virtual_machines.power_off(self.GROUP_NAME, self.VM_NAME)
        async_vm_stop.wait()

    def update_disk(self, increment=128):
        print("Expanding OS disk.")
        virtual_machine = self.compute_client.virtual_machines.get(self.GROUP_NAME,self.VM_NAME)
        os_disk_name = virtual_machine.storage_profile.os_disk.name
        os_disk = self.compute_client.disks.get(self.GROUP_NAME, os_disk_name)
        os_disk.disk_size_gb += increment
        async_disk_update = self.compute_client.disks.create_or_update(self.GROUP_NAME, os_disk.name, os_disk)
        async_disk_update.wait()

    def setup_docker(self):
        print("Loading toolkit.")
        run_command_parameters = {
          'command_id': 'RunShellScript', # For linux, don't change it
          'script': [
              'sudo apt install docker.io'
          ]
        }
        poller = self.compute_client.virtual_machines.run_command(
            self.GROUP_NAME,
            self.VM_NAME,
            run_command_parameters
        )
        result = poller.result()
        print(result.value[0].message)

    def check_status(self):
        pass

    def delete(self):
        print("Terminating.")
        async_vm_delete = self.compute_client.virtual_machines.delete(self.GROUP_NAME, self.VM_NAME)
        async_vm_delete.wait()

    def get_credentials(self):
        subscription_id = os.environ['AZURE_SUBSCRIPTION_ID']
        credentials = ServicePrincipalCredentials(
            client_id=os.environ['AZURE_CLIENT_ID'],
            secret=os.environ['AZURE_CLIENT_SECRET'],
            tenant=os.environ['AZURE_TENANT_ID']
        )
        return credentials, subscription_id

if __name__ == '__main__':
    new_vm = VirtualMachine()
    new_vm.deallocate()
    new_vm.update_disk()
    new_vm.start()
    new_vm.setup_docker()
    new_vm.delete()
