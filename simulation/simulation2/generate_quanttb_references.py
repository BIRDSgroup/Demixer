import os
import numpy as np
import random
import sys
import pandas as pd
import pdb

#mkdir ref_quanttb
def read_reference():
    with open('ref/CP003248.2.fasta', 'r') as refFile:
        next(refFile)
        refSeq = refFile.read().replace('\n', '')
    return list(refSeq)

print("reading SNPs")
snps=pd.read_csv('db/quanttb_snps.txt', delimiter=' ',header=None)
print("reading strain names")
lineages=pd.read_csv('db/quanttb_lineages.txt', delimiter=' ',header=None)
print("reading snp positions")
lineage_snppos=pd.read_csv('db/quanttb_snps.txt', delimiter=' ',header=None)
snp_matrix=np.load("db/snpmatrix.npy")

row,col=snp_matrix.shape

#print(row,col)
for i in range(0,2167):

    #lin1=lineages[i].values[0]
    
    
    lin_name=lineages[i].values[0]
    strain=f'>gi|{lin_name}|'
    

    x=((np.where(snp_matrix[i]==1)))[0]
    
    y=[int(lineage_snppos[v].values[0][:-1]) for v in x ]
    
    refSeq = read_reference()
    
    for j in range(0,len(y)):
        #thaprint(lin_name,j,y[j])
        refSeq[y[j]-1]=lineage_snppos[x[j]].values[0][-1]
    
    
    with open('ref_quanttb/'+str(lin_name)+'.fasta', 'w') as f:
    
        f.write(strain + '\n')
        str1="";
        str1=str1.join(refSeq)
        f.write(str1)
    
    


