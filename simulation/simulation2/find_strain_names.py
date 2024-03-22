import random
import sys
import pandas as pd
import warnings

#input arguments strain_number and db_name
seed = int(28082023)
random.seed(seed) 
strain_num=int(sys.argv[1]);

#lineage_names=pd.read_csv('TBprofiler_snpset_old_new.csv').columns
warnings.filterwarnings("ignore")
if(sys.argv[2]=="tbprof"): 
    lineage_names=pd.read_csv('../../scripts/db/tbprof_snpset.csv').columns
    lineage_snps=pd.read_csv('../../scripts/db/tbprof_unique.csv')
else:
    lineage_names=pd.read_csv('../../scripts/db/quanttb_snpset.csv').columns
    lineage_snps=pd.read_csv('../../scripts/db/quanttb_unique.csv')

for i in range(0,200):
    lin_id=[];
    strain_names="";
    for i in range(0,strain_num):
        
        flag=1;
        
        if(sys.argv[2]=="tbprof"): 
            while(flag==1):
                ind=random.randint(0,90)
                x=lineage_snps[lineage_names[ind]];
                y=[v for v in x if v == v]
                #print(len(x),len(y))
                if(len(y)>5 and lineage_names[ind]!='4.9' and lineage_names[ind]!='4'):
                    flag=0;
            strain_names=strain_names+'lineage'+lineage_names[ind]+'.fasta' + "\t"
        
        else:
            while(flag==1):
                ind=random.randint(0,2165)
                x=lineage_snps[lineage_names[ind]];
                y=[v for v in x if v == v]
                if(len(y)>5):
                    flag=0;
            strain_names=strain_names+lineage_names[ind]+'.fasta' + "\t"
        
    print(strain_names);


