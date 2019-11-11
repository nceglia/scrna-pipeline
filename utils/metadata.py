import pickle
import os
import collections
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
import pprint

SCOPES                = ['https://www.googleapis.com/auth/spreadsheets.readonly']
SAMPLE_SPREADSHEET_ID = '1plhIL1rH2IuQ8b_komjAUHKKrnYPNDyhvNNRsTv74u8'
SAMPLE_RANGE_NAME     = 'sample_metadata!A1:L'

pp = pprint.PrettyPrinter(indent=2)

def all_samples(index_by="nick_unique_id"):
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
        samples[sample[index_by]] = sample
    return samples

def lookup_by_igo(igo_id):
    samples = all_samples(index_by="unique_id")
    valid_sample = None
    if igo_id in samples:
        valid_sample = samples[igo_id]
    return valid_sample

def patient_samples(patient_id):
    samples = all_samples()
    valid_samples = []
    for unique_id, fields in samples.items():
        if fields["patient_id"] == patient_id:
            valid_samples.append(fields)
    return valid_samples

def single_sample(nick_unique_id):
    return samples()[nick_unique_id]

if __name__ == '__main__':
    all_samples = all_samples()
    pp.pprint(all_samples)
