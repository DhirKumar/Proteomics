#!/usr/bin/env python
# coding: utf-8

# In[49]:


import random
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
    

