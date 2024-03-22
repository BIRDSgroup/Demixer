import os
import numpy as np
import random
import sys
import pandas as pd


#mkdir ref_tbprof
def read_reference():
    with open('ref/tuberculosis.fna', 'r') as refFile:
        next(refFile)
        refSeq = refFile.read().replace('\n', '')
    return list(refSeq)

lineage_seedwords=pd.read_csv('../../scripts/db/tbprof_unique.csv')
#lineage_seedwords=pd.read_csv('TBprofiler_snpset_v3_new.csv')
#lineage_names=pd.read_csv('TBprofiler_snpset_old_new.csv').columns

lineage_names=pd.read_csv('../../scripts/db/tbprof_snpset.csv').columns

#lineage_names=pd.read_csv('TBprofiler_snpset_old_new.csv').columns

for i in range(0,91):

    lin1=lineage_seedwords.columns[i]
    lin_name='lineage'+lineage_names[i]
      
    strain=f'>gi|{lin_name}|'
    x=lineage_seedwords[lin1].values.flatten()
    
    for v in x:
        if v == v and v!='C4309264G,T' and v!='C786815G,T':
            print(v);
            print(int(v[1:-1]))

    y=[int(v[1:-1]) for v in x if v == v and v!='C4309264G,T' and v!='C786815G,T']

    
    
    refSeq = read_reference()
    
    for i in range(0,len(y)):
        if(x[i][-1]!=refSeq[y[i]]):
            refSeq[y[i]-1]=x[i][-1]
    
    
    with open('ref_tbprof/'+lin_name+'.fasta', 'w') as f:
    
        f.write(strain + '\n')
        str1="";
        str1=str1.join(refSeq)
        f.write(str1)
    
    


