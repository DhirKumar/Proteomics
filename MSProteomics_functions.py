#!/usr/bin/env python
# coding: utf-8

# In[33]:


import random
from pyteomics import mzid
import pandas as pd

def readFasta(fastapath):
    fastadict={}
    idx=''
    seq=''
    with open(fastapath,"r") as ff:
        for line in ff.readlines():
            if line.startswith('>'):
                if seq != '':
                    fastadict[idx]=seq
                idx=line.split(' ')[0].replace('>','')
                seq=''
            else:
                seq+=line.strip().replace(' ','')
    return(fastadict)
def createReverseDecoy(fastadict):
    decoydict={}
    for key in fastadict.keys():
        dkey='Decoy_'+key
        seq=fastadict[key][::-1]
        decoydict[dkey]=seq
    combineddict=fastadict | decoydict
    return(combineddict)

def createRandomizedDecoy(fastadict):
    decoydict={}
    for key in fastadict.keys():
        dkey='Decoy_'+key
        seqlist=list(fastadict[key])
        random.shuffle(seqlist)
        seq=''.join(seqlist)
        decoydict[dkey]=seq
    combineddict=fastadict | decoydict
    return(combineddict)

def extract_idents_mzid(mzidfile):
    records=[]
    with mzid.read(mzidfile) as reader:
        for item in reader:
            spectrum_id=item['spectrumID']
            rt=item['retention time']
            for sii in item['SpectrumIdentificationItem']:
                pep = sii.get('PeptideSequence')
                pass_thresh=sii.get('passThreshold')
                rnk=sii.get('rank')
                mods=[]
                mod_cc=''
                if sii.get('Modification'):
                    for mod in sii.get('Modification'):
                        #mods.append(mod.get('name'))
                        mods.append('_'.join([mod.get('name'),str(mod.get('location'))]))
                    mod_cc=';'.join(mods)
            records.append({'spec_id':spectrum_id,'RT':rt,'Peptide':pep,'rank':rnk,'mod':mod_cc,'pass_thresh':pass_thresh})
    return(records)

