#Code for generating the 300 read sequences for comparing Demixer with DL models 

import sys
from vcf_preprocess import*
import pickle
import pdb
import os
from Bio import SeqIO
import pandas as pd
import subprocess
import random
import csv

random.seed(42)
file_path = './db/NC_000962.3.fasta'

# Reading the FASTA file
seq_record = list(SeqIO.parse(file_path, "fasta"))[0]
ref_fasta=''.join(list(seq_record.seq))

path="../simulation/simulation2/ref_tbprof/"
l_path=[path+"lineage1.fasta",path+"lineage2.fasta",path+"lineage3.fasta",path+"lineage4.fasta"]



with open('finaloutput/Cryptic_rerun/Docs_col_name.pkl', 'rb') as f:
    setargs.Docs_col_name=pickle.load(f)
f.close();

setargs.Docs_test=np.zeros((300,148990),dtype="int16")
setargs.cols_name=dict(zip(setargs.Docs_col_name,range(len(setargs.Docs_col_name))))

j=0
values = {10, 25, 35, 40, 60}
lineage=[]

with open("Deep_train/random_lineage.fasta", 'w') as fasta_file:
    for l in l_path:
        seq_record = list(SeqIO.parse(l, "fasta"))[0]
        l1_fasta=''.join(list(seq_record.seq))
        for i in range(0,len(l1_fasta)):
            if(l1_fasta[i]!=ref_fasta[i]):
                col_id=ref_fasta[i]+str(i+1)+l1_fasta[i]
                if(setargs.cols_name.get(col_id) is not None):
                    ind=setargs.cols_name.get(col_id)
                    setargs.Docs_test[j,ind]=1
                    j+=1
                    random_value = random.choice(list(values))
                    start_position=i-random_value
                    substring = l1_fasta[start_position:start_position + 100]
                    name=">"+l.split("/")[-1][:-6]
                    fasta_file.write(f"{name}\n")
                    fasta_file.write(f"{substring}\n")
                    lineage.append((name[1:]))

with open('Deep_train/lineage.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(lineage)

        
out_dir = 'finaloutput/test'
try:
    os.mkdir(out_dir)
except OSError as error:
    print(error)      

filename = "finaloutput/test/Docs_test.dat"
fileobj = open(filename, mode='wb')
setargs.Docs_test.tofile(fileobj)
fileobj.close
        
with open('finaloutput/test/'+'testconfig.txt', 'w') as file:
    file.write(str(len(setargs.Docs_test)));

