#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 09:35:27 2024

@author: jlehle
"""
#%%
import os, sys, tarfile
os.chdir(os.environ['HOME']+'/OmniScrape/tmp')
#%%
def extract(tar_url, extract_path='.'):
    print(tar_url)
    tar = tarfile.open(tar_url, 'r')
    for item in tar:
        tar.extract(item, extract_path)
        if item.name.find(".tgz") != -1 or item.name.find(".tar") != -1:
            extract(item.name, "./" + item.name[:item.name.rfind('/')])
try:

    extract(sys.argv[1] + '.tgz')
    print('Done.')
except:
    name = os.path.basename(sys.argv[0])
    print(name[:name.rfind('.')], '<filename>')
#%%
extract('GSE252723_family.xml.tgz')

#%%
import requests
import re
import csv
import json
import time
from random import seed
from random import randint
#%%
from bs4 import BeautifulSoup 
# Reading the data inside the xml
# file to a variable under the name 
# data
#%%
with open('GSE252723_family.xml', 'r') as f:
    data = f.read()
    
#%%
# Passing the stored data inside
# the beautifulsoup parser, storing
# the returned object 
Bs_data = BeautifulSoup(data, "xml")
print(Bs_data)
#%%
import xml.etree.ElementTree as ET
xmlTree = ET.parse('GSE252723_family.xml')
print(xmlTree)
tags = {elem.tag for elem in xmlTree.iter()}
tags
#%%
# Finding all instances of tag 
# `unique`
b_unique = Bs_data.find_all('Data-Table')
b_unique
#%%
headers = {
  'authority': 'screen-beta-api.wenglab.org',
  'accept': 'application/json',
  'accept-language': 'en-US,en;q=0.9',
  'content-type': 'application/json',
  'origin': 'https://screen.wenglab.org',
  'referer': 'https://screen.wenglab.org/',
  'sec-ch-ua': '\".Not/A)Brand\";v=\"99\", \"Google Chrome\";v=\"103\", \"Chromium\";v=\"103\"',
  'sec-ch-ua-mobile': '?1',
  'sec-ch-ua-platform': '\"Android\"',
  'sec-fetch-dest': 'empty',
  'sec-fetch-mode': 'cors',
  'sec-fetch-site': 'same-site',
  'user-agent': 'Mozilla/5.0 (Linux; Android 6.0; Nexus 5 Build/MRA58N) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/103.0.0.0 Mobile Safari/537.36'
}
#%%
url = 'https://www.ncbi.nlm.nih.gov/gds/?term=200264681'
#%%
url_output = requests.get(url)
url_output.headers['Content-Type']
url_output.text
geo_json = url_output.json()

pd.readhttps://www.ncbi.nlm.nih.gov/gds/?term=200264681
type(X)
X
X.loc[0,0].startswith('1.')


import requests
fname = 'guppy-0.1.10.tar.gz'
url = 'https://pypi.python.org/packages/source/g/guppy/' + fname
r = requests.get(url)
open(fname , 'wb').write(r.content)

#%%
#######
# THIS IS WHERE THINGs STARTED WORKING #
#######
#%%
#%%
import pandas as pd
X = pd.read_fwf('/master/jlehle/OmniScrape/tmp/gds_result.txt', header=None)
ftp_links = X[X[0].str.startswith('FTP download')]
ftp_links = ftp_links.reset_index(drop=True)
ftp_links
#%%
ftp_list = []
# loop through the rows using iterrows()
for index, row in ftp_links.iterrows():
    try:
        tmp = row.to_string()
        tmp = tmp.split(') ', -1)[1]
        ftp_list.append(tmp)
    except:
        pass

len(ftp_list)

#%%
import GEOparse
#%%
gse = GEOparse.get_GEO(filepath="GSE202695_family.soft.gz")
gse = GEOparse.get_GEO(filepath="GSE161529_family.soft.gz")
gse.gpls
gse.gsms
gse.gsms['GSM4909314'].show_metadata()
gse.gsms['GSM6129415'].get_metadata_attribute('Sample_title')
gse.gsms['GSM4909314'].download_SRA('jlehle@txbiomed.org', directory=GSE_UID, nproc=1,)
gse.gsms['GSM4909314'].download_supplementary_files()

pd.set_option('display.max_columns', None)
print(gse.phenotype_data)
print(gse.phenotype_data['geo_accession'])

for accession in gse.phenotype_data['geo_accession'];:
    gse.gsms[accession].download_SRA('jlehle@txbiomed.org', directory=gse.phenotype_data['series_id'], nproc=1)
print(gse.phenotype_data['series_id'])
