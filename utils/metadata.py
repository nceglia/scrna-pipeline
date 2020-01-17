import pickle
import os
import collections
import urllib
import requests
import glob

from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request

import pprint

API_KEY     = os.environ["ELAB_API_KEY"]
SAMPLES_URL = "https://msk.elabinventory.com/api/v1/samples/"
SAMPLE_URL  = "https://msk.elabinventory.com/api/v1/samples/{sampleid}/meta"

SCOPES                = ['https://www.googleapis.com/auth/spreadsheets.readonly']
SAMPLE_SPREADSHEET_ID = '1plhIL1rH2IuQ8b_komjAUHKKrnYPNDyhvNNRsTv74u8'
SAMPLE_RANGE_NAME     = 'sample_metadata!A1:L'

pp = pprint.PrettyPrinter(indent=2)

def resolve_igo(unique_id):
    validate="/work/shah/data/scrnaseq/*/"
    valid_directories = glob.glob(validate)
    for valid_directory in valid_directories:
        valid_sample = valid_directory.split("/")[-2]
        if "_".join(valid_sample.split("_")[:-1]) == unique_id:
            return valid_sample
    return None

def all_samples_elab(index_by="nick_unique_id"):
    index_samples = dict()
    headers = {'Authorization': API_KEY, "Host": "msk.elabinventory.com"}
    response = requests.get(SAMPLES_URL, headers=headers)
    samples = response.json()
    for sample in samples["data"]:
        sampleid = sample["sampleID"]
        patient = sample["name"]
        if patient == "SPECTRUM-OV-000": continue
        response = requests.get(SAMPLE_URL.format(sampleid=sampleid), headers=headers)
        sampleinfo = response.json()
        metadata = dict([(meta["key"],meta["value"]) for meta in sampleinfo['data'] if "value" in meta])
        if "scRNA IGO ID" in metadata:
            igo_project = metadata["scRNA IGO ID"]
        elif "IGO ID" in metadata:
            continue
        if "scRNA REX Sample ID" in metadata:
            rex_prefix  = metadata["scRNA REX Sample ID"]
        elif "REX Sample ID" in metadata:
            continue
        if igo_project != "" and rex_prefix != "":
            igo = "Sample_{rex_prefix}_IGO_09443_{igo_project}".format(rex_prefix=rex_prefix, igo_project=igo_project)
            igo = resolve_igo(igo)
            if not igo:
                continue
            subsite = metadata["Site Details"]
            site = metadata["Specimen Site"]
            if subsite.strip() == '':
                subsite = site
            sort = metadata["Submitted Populations"].replace("+","P").replace("-","N")
            surgery = "S{}".format(int(metadata["Surgery #"])+1)
            internal = "{patient}_{surgery}_{sort}_{site}".format(patient=patient, sort=sort, surgery=surgery, site="_".join(subsite.split())).upper()
            if index_by.lower() == 'unique_id':
                primary_key = igo
            else:
                primary_key = internal
            relevant_metadata = {"site": site.upper(), "surgery": surgery, "sort": sort, "unique_id": igo, "nick_unique_id": internal, "patient_id": patient, "subsite": "_".join(subsite.upper().split()), "origin": "elab"}
            index_samples[primary_key] = relevant_metadata
    return index_samples

def all_samples_gs(index_by="nick_unique_id"):
    creds = None
    if os.path.exists('token.pickle'):
        with open('token.pickle', 'rb') as token:
            creds = pickle.load(token)
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file('credentials.json', SCOPES)
            creds = flow.run_local_server(port = 0)
        with open('token.pickle', 'wb') as token:
            pickle.dump(creds, token)

    service = build('sheets', 'v4', credentials = creds)
    sheet = service.spreadsheets()
    result = sheet.values().get(spreadsheetId = SAMPLE_SPREADSHEET_ID,
                                range = SAMPLE_RANGE_NAME).execute()
    values = result.get('values', [])

    samples = collections.defaultdict(dict)
    header = values.pop(0)
    for row in values:
        sample = dict(zip(header,row))
        sample["origin"] = "google_sheet"
        samples[sample[index_by]] = sample
    return samples

def all_samples(index_by='nick_unique_id'):
    samples = all_samples_elab(index_by=index_by)
    gs_samples = all_samples_gs(index_by=index_by)
    for sample in gs_samples:
        if sample not in samples:
            samples[sample] = gs_samples[sample]
    return samples

def lookup_by_igo(igo_id):
    samples = all_samples(index_by="unique_id")
    if igo_id in samples:
        valid_sample = samples[igo_id]
    elif "_".join(igo_id.split("_")[:-1]) in samples:
        valid_sample = samples["_".join(igo_id.split("_")[:-1])]
    return valid_sample

def patient_samples(patient_id):
    samples = all_samples(index_by='nick_unique_id')
    valid_samples = []
    for unique_id, fields in samples.items():
        if fields["patient_id"] == patient_id:
            valid_samples.append(fields)
    return valid_samples

def single_sample(nick_unique_id):
    samples = all_samples(index_by='nick_unique_id')
    sample = samples.get(nick_unique_id,None)
    return sample

if __name__ == '__main__':
    samples = all_samples()
    pp.pprint(samples)
    origin = collections.defaultdict(int)
    patients = collections.defaultdict(int)
    for sample, metadata in samples.items():
        origin[metadata["origin"]] += 1
        patients[metadata["patient_id"]] += 1
    for source, count in origin.items():
        print(source, count)
    print(len(patients.items()))
    for patient in sorted(list(patients.keys())):
        print(patient,patients[patient])
    valid_samples = patient_samples("SPECTRUM-OV-054")
    pp.pprint(valid_samples)
    valid_sample = lookup_by_igo("Sample_054RA_CD45P_IGO_09443_AC_2")
    print(valid_sample)
