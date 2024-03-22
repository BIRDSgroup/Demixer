from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd;
import allel
import numpy as np
import sys
from operator import add 
import os
current_directory = os.getcwd()
two_levels_up = os.path.abspath(os.path.join(current_directory, "../../scripts"))
sys.path.append(two_levels_up)

import setargs

ppe=pd.read_excel('../../scripts/db/12864_2020_6486_MOESM4_ESM.xlsx',index_col=0)

#print(ppe)
ppe_start=ppe['start'].values;
ppe_end=ppe['end'].values;
#print(genotype[((genotype['POS']>=ppe_start[3])&(genotype['POS']<=ppe_end[3])==True)])



def process_row(i,ref,length,arg=2):
    ad=list();
    
    

    
    if(arg==1):
        if(ref==0):
            ad.extend(setargs.callset['calldata/RO'][i]);
            ad=[0 if x==-1 else x for x in ad]
        else:
            ad.extend(setargs.callset['calldata/AO'][i][:,ref-1]);
            ad=[0 if x==-1 else x for x in ad]
        #print(ad)
    elif(arg==2):
    
    
    
        for j in range(0,len(setargs.callset['samples'])): #no. of samples
            
            if(setargs.callset['calldata/AD'][i][j][ref]>0):
                ad.append(setargs.callset['calldata/AD'][i][j][ref])
            else:
                ad.append(0)
    else:
        for j in range(0,len(setargs.callset['samples'])): #no. of samples
            
            if(setargs.callset['calldata/COV'][i][j][ref]>0):
                ad.append(setargs.callset['calldata/COV'][i][j][ref])
            else:
                ad.append(0)
        
    
    return ad[0:length];

def cal_genotype(filename):
    start=[];
    end=[];
    for  line in open("../../scripts/Mycobacterium_tuberculosis_annotation.gff3"):
        if not line.lstrip().startswith('#'):
            line_arr = line.split("\t")
            if(line_arr[2]=="CDS"):
                #print(line)
                start.append(int(line_arr[3]));
                end.append(int(line_arr[4]));
    genotype=allel.vcf_to_dataframe(filename,fields=['POS','REF','ALT','is_snp']);  
    
    print('hello');
    for i in range(0,len(start)-1):  
        genotype=genotype[((genotype['POS']>end[i])&(genotype['POS']<start[i+1])==False)]
        #print(len(genotype.index))

    for i in range(0,len(ppe_start)-1):  
        genotype=genotype[((genotype['POS']>=ppe_start[i])&(genotype['POS']<=ppe_end[i])==False)]
        #print(len(genotype.index))
    #print(genotype)
    return genotype;





def read_data(i,id_num,genotype,length,read_countdata):
   
    #print(genotype['POS'],setargs.callset['variants/POS'])
   
    pos_index = (np.where(setargs.callset['variants/POS']==genotype['POS'][i]))[0][0]
    #print(pos_index,i)

    if(pos_index>=0):
        
        if(not pd.isnull(genotype[('ALT_'+str(1))].iat[i])):
            col_name=genotype['REF'].iat[i]+str(genotype['POS'].iat[i])  + str(genotype[('ALT_'+str(1))].iat[i])
            ad=process_row(pos_index,1,length,1);
            
 
            if(1):
        
                col_name=genotype['REF'].iat[i]+str(genotype['POS'].iat[i]) +genotype['REF'].iat[i]

                if(col_name not in setargs.Docs_col_name):     
                        ad=process_row(pos_index,0,length,1);

                        setargs.Docs_col_name.append(col_name);
                      
                        for kk in range(length):
                            if(ad[kk]>5 ):
                                
                                setargs.Docs_1[kk].append(ad[kk])
 
                            else:
                                setargs.Docs_1[kk].append(0)

                        id_num=id_num+1;

                        

                        for j in range(1,4):
                            #print(j,id_num)
                            if(not pd.isnull(genotype[('ALT_'+str(j))].iat[i])):
                                col_name=genotype['REF'].iat[i]+str(genotype['POS'].iat[i])  + str(genotype[('ALT_'+str(j))].iat[i])
                                
                                ad=process_row(pos_index,j,length,1);
 
                                setargs.Docs_col_name.append(col_name);

                                for kk in range(length):
                                    if(ad[kk]>5):
                                        #setargs.Docs[kk].extend(np.full(ad[kk],id_num,dtype=int))
                                        setargs.Docs_1[kk].append(ad[kk])
                                    else:
                                        setargs.Docs_1[kk].append(0)

                                id_num=id_num+1;
                            else:
                                    col_name=col_name=genotype['REF'].iat[i]+str(genotype['POS'].iat[i]) + 'ALT'+str(j)
                                    setargs.Docs_col_name.append(col_name)
                                    for kk in range(length):
                                        setargs.Docs_1[kk].append(0)

                                    id_num=id_num+1;
        
        return id_num;







