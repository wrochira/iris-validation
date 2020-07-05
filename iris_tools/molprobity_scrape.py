"""
Copyright 2020 William Rochira at York Structural Biology Laboratory

- I wrote this before I discovered the MMTBX Python Molprobity module
- It's probably rude to scrape the Duke server for every single PDB model
- So don't use it
- But here it is
"""

import os
import time
import requests
from random import shuffle

try:
    from BeautifulSoup import BeautifulSoup
except ImportError:
    from bs4 import BeautifulSoup

from _defs import MOLPROBITY_OUTPUT_DIR
from common import setup, get_available_pdb_ids


MOLPROBITY_URL = 'http://molprobity.biochem.duke.edu/'
MOLPROBITY_VERSION = 'MolProbity20200513_4.5.1'


def new_molprob_session():
    response = requests.get(MOLPROBITY_URL)
    body = BeautifulSoup(response.text, 'html.parser').body
    toplink = body.a['href']
    msid = toplink.split('MolProbSID=')[1].split('&')[0]
    return msid

def post_pdb(msid, pdb_id):
    post_data = { 'MolProbSID' : msid,
                  'eventID' : '14',
                  'pdbCode' : pdb_id,
                  'fetchType' : 'pdb',
                  'cmd' : 'Fetch >',
                  'uploadFile' : '',
                  'uploadType' : 'pdb' }
    response = requests.post(MOLPROBITY_URL + 'index.php', data=post_data)

def job_finished(msid):
    params = { 'MolProbSID' : msid }
    response = requests.get(MOLPROBITY_URL + 'index.php', params=params)
    head = BeautifulSoup(response.text, 'html.parser').head
    if head is None:
        return False
    if 'Job is finished' in head.title.text:
        return True
    return False

def acknowledge_submission(msid):
    params = { 'MolProbSID' : msid,
               'eventID' : 23 }
    requests.get(MOLPROBITY_URL + 'index.php', params=params)
    post_data = { 'MolProbSID' : msid,
                  'eventID' : '24',
                  'cmd' : 'Continue >' }
    response = requests.post(MOLPROBITY_URL + 'index.php', data=post_data)
    params = { 'MolProbSID' : msid,
               'eventID' : 44 }
    requests.get(MOLPROBITY_URL + 'index.php', params=params)

def initiate_analysis(msid, pdb_id):
    post_data = { 'MolProbSID' : msid,
                  'eventID' : '61',
                  'modelID' : pdb_id,
                  'doKinemage' : '1',
                  'kinGeom' : '1',
                  'kinRama' : '1',
                  'kinRota' : '1',
                  'kinCBdev' : '1',
                  'kinOmega' : '1',
                  'doCharts' : '1',
                  'chartGeom' : '1',
                  'chartRama' : '1',
                  'chartRota' : '1',
                  'chartCBdev' : '1',
                  'chartOmega' : '1',
                  'chartMulti' : '1',
                  'chartNotJustOut' : '1',
                  'chartAltloc' : '1',
                  'cmd' : 'Run programs to perform these analyses >' }
    response = requests.post(MOLPROBITY_URL + 'index.php', data=post_data)

def acknowledge_analysis(msid):
    params = { 'MolProbSID' : msid,
               'eventID' : 62 }
    requests.get(MOLPROBITY_URL + 'index.php', params=params)

def get_multi_table(msid, pdb_id):
    table_values = [ ]
    #url = MOLPROBITY_URL + 'viewtable.php?MolProbSID=' + msid + '&file=/home/htdocs/molprobity/public_html/data/' + msid + '/raw_data/' + pdb_id + '-multi.table'
    url = MOLPROBITY_URL + 'viewtable.php?MolProbSID=' + msid + '&file=/mental-data/rlab/MolProbities/' + MOLPROBITY_VERSION + '/molprobity/public_html/data/' + msid + '/raw_data/' + pdb_id + '-multi.table'
    response = requests.get(url)
    html = BeautifulSoup(response.text, 'html.parser')
    if html is None:
        print('Bad result, check URL definition.')
        return None
    try:
        table = html.find_all('table')[-1]
        relevant_rows = table.find_all(lambda tag: tag.name=='tr' and tag['bgcolor'] in ('#ffffff', '#f0f0f0'))
    except:
        print('Bad table')
        return None
    for i, row in enumerate(relevant_rows):
        cells = row.find_all('td')
        try:
            chainres, _, code, b_high, rama_info, rota_info, cb_dev, _, _, _ = [ cell.text.strip() for cell in cells ]
        except:
            print('Unexpected formatting on row', i+1)
            return None
        chainres = chainres.replace(' ', '')
        chain = ''.join([ x for x in chainres if not x.isdigit() ])
        resnum = chainres[len(chain):]
        b_high = float(b_high)
        phi = psi = rama_probability = None
        if rama_info != '-':
            phi, psi = [ float(x) for x in rama_info.split(' / ')[1].split(',') ]
            rama_probability = float(rama_info.split('(')[1].split('%')[0].replace(')', ''))
        chis = [ ]
        rotamer_name = sc_probability = None
        if rota_info != '-':
            chis = [ float(x) for x in rota_info.split('chi angles: ')[1].split(',') ]
            rotamer_name = rota_info.split('chi angles: ')[0].split(' ')[-1] if 'OUTLIER' not in rota_info else None
            sc_probability = float(rota_info.split('(')[1].split('%')[0])
        chis = (chis + [ None ] * 4)[:4]
        table_values.append((chain, resnum, code, b_high, phi, psi, rama_probability, *chis, rotamer_name, sc_probability))
    return table_values



if __name__ == '__main__':
    setup()
    pdb_ids = get_available_pdb_ids()
    for pdb_id in pdb_ids:
        pdb_id = pdb_id.lower()
        print()
        print('PDB ID:', pdb_id)
        outfile_name = MOLPROBITY_OUTPUT_DIR + pdb_id + '.csv'
        if os.path.exists(outfile_name):
            print('Skipping', pdb_id)
        print('Starting new Molprobity session...')
        msid = new_molprob_session()
        print('MSID:', msid)
        print('Submitting PDB code', pdb_id, 'for analysis...')
        post_pdb(msid, pdb_id)
        print('Waiting for submission to complete...')
        while not job_finished(msid):
            time.sleep(1)
        print('Acknowledging receipt of submission...')
        acknowledge_submission(msid)
        print('Initiating AAC and geometry analysis...')
        initiate_analysis(msid, pdb_id)
        print('Waiting for analysis to complete...')
        while not job_finished(msid):
            time.sleep(1)
        print('Acknowledging completion of analysis...')
        acknowledge_analysis(msid)
        print('Getting multicritereon table...')
        table_values = get_multi_table(msid, pdb_id)
        if table_values is None:
            continue
        print('Writing to output file...')
        with open(outfile_name, 'w') as outfile:
            outfile.write('Chain ID, Sequence Number,Residue Code,Max B,Phi,Psi,Ramachandram Probability,Chi1,Chi2,Chi3,Chi4,Rotamer Name,Rotamer Probability\n')
            for line in table_values:
                outfile.write(','.join(str(x) for x in line) + '\n')
        print('Done.')
    print()
