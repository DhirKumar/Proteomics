#!/usr/bin/env python
# coding: utf-8
import pandas as pd
fasta_file='uniprotkb_proteome_UP000005640_AND_revi_2024_05_24.fasta'
motifs=['RKH','HRK','HKR','RHK','KRH','KHR']

def read_fasta(fasta_file):
    count=0
    seq_dict={}
    seq=''
    idx=''
    with open(fasta_file,'r') as ff:
        for line in ff.readlines():
            line=line.strip()
            if(line.startswith('>')) and seq != '':
                seq_dict[idx]=seq
                idx= line.split(' ')[0].replace('>','')
                seq=''
            else:
                seq+=line.strip()
    return(seq_dict)

def motif_location(seq,motif):
    return(seq.find(motif))

all_seq_dict=read_fasta(fasta_file)

motif_mapped=dict()
for idx in all_seq_dict.keys():
    for motif in motifs:
        loc=motif_location(all_seq_dict[idx],motif)
        if loc >= 0:
            #print(loc)
            if motif_mapped.get('idx','') =='':
                motif_mapped[idx]={motif:loc}
                
df=pd.DataFrame(motif_mapped)
df.head()

   




