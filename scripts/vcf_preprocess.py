from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd;
import allel
import numpy as np
import setargs;
import sys
import pdb

ppe=pd.read_excel('db/12864_2020_6486_MOESM4_ESM.xlsx',index_col=0)

#print(ppe)
ppe_start=ppe['start'].values;
ppe_end=ppe['end'].values;
#print(genotype[((genotype['POS']>=ppe_start[3])&(genotype['POS']<=ppe_end[3])==True)])

variants=pd.read_csv('db/'+sys.argv[3]+'_snpset.csv')        #for invitro
 
#variants=pd.DataFrame()   #included for ablation study - to be removed

seeds=variants.stack().values    

def process_row(i,ref,length,arg=2):
    ad=list();
    if(arg==1):
        if(ref==0):
            ad.extend(setargs.callset['calldata/RO'][i]);
            ad=[0 if x==-1 else x for x in ad]
        else:
            ad.extend(setargs.callset['calldata/AO'][i][:,ref-1]);
            ad=[0 if x==-1 else x for x in ad]
     
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
    for  line in open("db/Mycobacterium_tuberculosis_annotation.gff3"):
        if not line.lstrip().startswith('#'):
            line_arr = line.split("\t")
            if(line_arr[2]=="CDS"):
                #print(line)
                start.append(int(line_arr[3]));
                end.append(int(line_arr[4]));
    genotype=allel.vcf_to_dataframe(filename,fields=['POS','REF','ALT','is_snp']);  
    
    
    for i in range(0,len(start)-1):  
        genotype=genotype[((genotype['POS']>end[i])&(genotype['POS']<start[i+1])==False)]
        #print(len(genotype.index))

    for i in range(0,len(ppe_start)-1):  
        genotype=genotype[((genotype['POS']>=ppe_start[i])&(genotype['POS']<=ppe_end[i])==False)]
        #print(len(genotype.index))
    #print(genotype)
    return genotype;


def read_data(i,genotype,length,read_countdata,weight,min_read=5):
   
    pos_index = (np.where(setargs.callset['variants/POS']==genotype['POS'][i]))[0][0]

    if(pos_index>=0):
        
        if(not pd.isnull(genotype[('ALT_'+str(1))].iat[i])):
            col_name1=genotype['REF'].iat[i]+str(genotype['POS'].iat[i])  + str(genotype[('ALT_'+str(1))].iat[i])
            col_name2=genotype['REF'].iat[i]+str(genotype['POS'].iat[i])  + str(genotype[('ALT_'+str(2))].iat[i])
            col_name3=genotype['REF'].iat[i]+str(genotype['POS'].iat[i])  + str(genotype[('ALT_'+str(3))].iat[i])
       
            col_name=genotype['REF'].iat[i]+str(genotype['POS'].iat[i]) +genotype['REF'].iat[i]
            
            if(1):   
                     
            #if((col_name in seeds) or (col_name1 in seeds) or (col_name2 in seeds) or (col_name3 in seeds)):
                             
                        ad=process_row(pos_index,0,length,int(sys.argv[4]));
                        
                        
                       
                        if(pd.isnull(genotype[('ALT_'+str(2))].iat[i])):
                            setargs.Docs_col_name1.append(col_name);
                            
                            if(col_name in seeds):
                                setargs.idf1.append(weight)
                            else:
                                setargs.idf1.append(0)
                  
                                        
                            for kk in range(length):
                                if(ad[kk]>min_read ):
                                
                                    setargs.Docs_1[kk].append(ad[kk])
                                
                                else:
                                    setargs.Docs_1[kk].append(0)
                            
                            if(not pd.isnull(genotype[('ALT_'+str(1))].iat[i])):
                                    col_name=genotype['REF'].iat[i]+str(genotype['POS'].iat[i])  + str(genotype[('ALT_'+str(1))].iat[i])

                                    ad=process_row(pos_index,1,length,int(sys.argv[4]));
                             
                                    setargs.Docs_col_name1.append(col_name);
                                    if(col_name in seeds):
                                        setargs.idf1.append(weight)
                                    else:
                                        setargs.idf1.append(0)
                       

                                    for kk in range(length):
                                        if(ad[kk]>min_read):
                                            
                                            setargs.Docs_1[kk].append(ad[kk])
                                        else:
                                            setargs.Docs_1[kk].append(0)

                            
                        elif(not pd.isnull(genotype[('ALT_'+str(2))].iat[i]) and pd.isnull(genotype[('ALT_'+str(3))].iat[i])):
                            
                            setargs.Docs_col_name2.append(col_name);
                            if(col_name in seeds):
                                setargs.idf2.append(weight)
                            else:
                                setargs.idf2.append(0)
                                       
                            for kk in range(length):
                                if(ad[kk]>min_read ):
                                
                                    setargs.Docs_2[kk].append(ad[kk])
                                
                                else:
                                    setargs.Docs_2[kk].append(0)
                            for j in range(1,3):
                                #print(j,id_num)
                                if(not pd.isnull(genotype[('ALT_'+str(j))].iat[i])):
                                    col_name=genotype['REF'].iat[i]+str(genotype['POS'].iat[i])  + str(genotype[('ALT_'+str(j))].iat[i])

                                    ad=process_row(pos_index,j,length,int(sys.argv[4]));
                                    #ad=process_row(pos_index,j,length,2)
                                    #print(col_name,ad);
                                    #read_countdata[col_name]=ad;
                                    setargs.Docs_col_name2.append(col_name);
                                    if(col_name in seeds):
                                        setargs.idf2.append(weight)
                                    else:
                                        setargs.idf2.append(0)


                                    for kk in range(length):
                                        if(ad[kk]>min_read):
                                            #setargs.Docs[kk].extend(np.full(ad[kk],id_num,dtype=int))
                                            setargs.Docs_2[kk].append(ad[kk])
                                        else:
                                            setargs.Docs_2[kk].append(0)

                            
                            
                        elif(not pd.isnull(genotype[('ALT_'+str(3))].iat[i])):
                            
                            setargs.Docs_col_name3.append(col_name);
                            
                            if(col_name in seeds):
                                setargs.idf3.append(weight)
                            else:
                                setargs.idf3.append(0)
   
            
                            for kk in range(length):
                                if(ad[kk]>min_read ):
                                
                                    setargs.Docs_3[kk].append(ad[kk])
                                
                                else:
                                    setargs.Docs_3[kk].append(0)

                            
                            for j in range(1,4):
                                #print(j,id_num)
                                if(not pd.isnull(genotype[('ALT_'+str(j))].iat[i])):
                                    col_name=genotype['REF'].iat[i]+str(genotype['POS'].iat[i])  + str(genotype[('ALT_'+str(j))].iat[i])

                                    ad=process_row(pos_index,j,length,int(sys.argv[4]));
                                    #ad=process_row(pos_index,j,length,2)
                                    #print(col_name,ad);
                                    #read_countdata[col_name]=ad;
                                    setargs.Docs_col_name3.append(col_name);

                                    if(col_name in seeds):
                                        setargs.idf3.append(weight)
                                    else:
                                        setargs.idf3.append(0)
                                    
                                    for kk in range(length):
                                        if(ad[kk]>min_read):
                                            
                                            setargs.Docs_3[kk].append(ad[kk])
                                        else:
                                            setargs.Docs_3[kk].append(0)
                            
def read_data_test(i,genotype,length,read_countdata,weight,min_read=5):
   
    pos_index = (np.where(setargs.callset['variants/POS']==genotype['POS'][i]))[0][0]

    if(pos_index>=0):
        
        if(not pd.isnull(genotype[('ALT_'+str(1))].iat[i])):
            
            col_name=genotype['REF'].iat[i]+str(genotype['POS'].iat[i]) +genotype['REF'].iat[i]
            ad=process_row(pos_index,0,length,int(sys.argv[4]));

            
            if(col_name in setargs.cols_name):

                ind=setargs.cols_name[col_name]

                
                if(pd.isnull(genotype[('ALT_'+str(2))].iat[i])):
                    for kk in range(length):
                        if(ad[kk]>min_read ):

                            setargs.Docs_test[kk][ind]=ad[kk]


                if(not pd.isnull(genotype[('ALT_'+str(1))].iat[i])):
                        col_name=genotype['REF'].iat[i]+str(genotype['POS'].iat[i])  + str(genotype[('ALT_'+str(1))].iat[i])
                        if(col_name in setargs.cols_name):
                            ind=setargs.cols_name[col_name]

                            ad=process_row(pos_index,1,length,int(sys.argv[4]));

                            for kk in range(length):
                                if(ad[kk]>min_read):

                                    setargs.Docs_test[kk][ind]=ad[kk]
                            
            elif(not pd.isnull(genotype[('ALT_'+str(2))].iat[i]) and pd.isnull(genotype[('ALT_'+str(3))].iat[i])):
                if(col_name in setargs.cols_name):

                    ind=setargs.cols_name[col_name]

                    for kk in range(length):
                        if(ad[kk]>min_read ):

                            setargs.Docs_test[kk][ind]=ad[kk]

                for j in range(1,3):
                    #print(j,id_num)
                    if(not pd.isnull(genotype[('ALT_'+str(j))].iat[i])):
                        col_name=genotype['REF'].iat[i]+str(genotype['POS'].iat[i])  + str(genotype[('ALT_'+str(j))].iat[i])
                        if(col_name in setargs.cols_name):
                            ind=setargs.cols_name[col_name]

                            ad=process_row(pos_index,j,length,int(sys.argv[4]));


                            for kk in range(length):
                                if(ad[kk]>min_read):
                                    #setargs.Docs[kk].extend(np.full(ad[kk],id_num,dtype=int))
                                    setargs.Docs_2[kk][ind]=ad[kk]                                          

            elif(not pd.isnull(genotype[('ALT_'+str(3))].iat[i])):

                    if(col_name in setargs.cols_name):

                        ind=setargs.cols_name[col_name]

                        for kk in range(length):
                            if(ad[kk]>min_read ):

                                setargs.Docs_test[kk][ind]=ad[kk]

                    for j in range(1,4):
                        #print(j,id_num)
                        if(not pd.isnull(genotype[('ALT_'+str(j))].iat[i])):
                            col_name=genotype['REF'].iat[i]+str(genotype['POS'].iat[i])  + str(genotype[('ALT_'+str(j))].iat[i])
                            
                            if(col_name in setargs.cols_name):
                                ind=setargs.cols_name[col_name]


                                ad=process_row(pos_index,j,length,int(sys.argv[4]));
   
                                for kk in range(length):
                                    if(ad[kk]>min_read):
                                        setargs.Docs_test[kk][ind]=ad[kk]










