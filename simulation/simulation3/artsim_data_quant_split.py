import re,sys
import os
import math
from operator import add 
import operator
import threading
import time
import numpy
from pandas import ExcelWriter
from pandas import ExcelFile
import scipy.sparse as sps
import pickle
import pdb;
import random
import pandas as pd
from preprocess import*
import numpy as np;
import multiprocessing
import matplotlib.pyplot as plt
from plotnine import *
from scipy import stats    
from sklearn.metrics import mean_absolute_error as mae
from scipy.stats import entropy
import time
current_directory = os.getcwd()
two_levels_up = os.path.abspath(os.path.join(current_directory, "../simulation1"))
sys.path.append(two_levels_up)

from init import*
from nmf import*

current_directory = os.getcwd()
two_levels_up = os.path.abspath(os.path.join(current_directory, "../../scripts"))
sys.path.append(two_levels_up)

from setargs import *
from init import *

rel_error=np.zeros((3,13))


current_directory = os.getcwd()
out_dir = os.path.join(current_directory, r'Demixer_cpp/input')
try:
    os.makedirs(out_dir)

except OSError as error:

    print(error)

out_dir = os.path.join(current_directory, r'Demixer_cpp/output')
try:
    os.makedirs(out_dir)
except OSError as error:
    print(error)

    
def var_initialize():
    
    setargs.z_m_n=[];
    #setargs.strain=sys.argv[3];    #output_filename
    setargs.n_m_z = numpy.zeros((len(setargs.Docs), setargs.K),dtype='int32')     # word count of each document and topic
    setargs.n_z_t = numpy.zeros((setargs.K, setargs.V),dtype='int32') # word count of each topic and vocabulary
    setargs.n_z = numpy.zeros(setargs.K)        # word count of each topic
    setargs.n_m_t = numpy.zeros((len(setargs.Docs), setargs.V,setargs.K),dtype='int32')
    setargs.n_r=numpy.zeros((len(setargs.Docs), int(setargs.V/4)+1))
    setargs.n_r_max=numpy.full((len(setargs.Docs), int(setargs.V/4)+1),-1)
    setargs.doc_term= numpy.zeros((len(setargs.Docs), setargs.V))

def mutations(n_z_t,est_theta):
    true_theta=pd.read_csv('artsim/'+str(no_strain)+'strain_true_theta.csv')
    
    
    for i in range(0,setargs.K):
        true_theta[str(i)]=true_theta[(str(i))]/100
        
    n_z_t=n_z_t+0.01
    phi1=np.zeros((setargs.K,setargs.V))
    for i in range(int(setargs.V/4)):
        n_z_t_sub=n_z_t[:,i*4:(i+1)*4]/n_z_t[:,i*4:(i+1)*4].sum(axis=1)[:,None]
        phi1[:,i*4:(i+1)*4]=numpy.nan_to_num(n_z_t_sub, nan=0)
    
    try:
        subset=pd.read_csv('all_strains.csv')
        row,col=subset.shape;
    except pd.errors.EmptyDataError:
        row=0
        col=0
    r,c=n_z_t.shape
    print(row,col)
    true_n_z_t = np.empty((col,c))
    true_n_z_t[:] = np.NaN
    l=len(subset.columns.values)
    
    for m in range(0,l):
            variant=list(subset.iloc[:,m])
            for i in range(0,len(variant)):
                    if(not math.isnan(variant[i])):
                        loc=int(variant[i]/4) 
                        #print(loc)
                        
                        true_n_z_t[m,loc*4:(loc+1)*4]=[0,1,0,0]
                                        
    diverg=["" for x in range(setargs.K)]
    
    
    for j in range(0,setargs.K):
            #print('j',j)
            divergence=np.zeros(l);
            for m in range(0,l):
                #print('m',m)
                variant=list(subset.iloc[:,m])
                count=0;
                for i in range(0,len(variant)):
                    if(not math.isnan(variant[i])):
                        #print(i,variant[i])
                        loc=int(variant[i]/4) 
                        divergence[m]=divergence[m]+entropy(true_n_z_t[m][loc*4:(loc+1)*4], phi1[j][loc*4:(loc+1)*4])
                        
                        count+=1
                divergence[m]=divergence[m]/count
            
            indice=np.where(divergence==np.nanmin(divergence))[0]
            print(np.nanmin(divergence),divergence)
            if(np.nanmin(divergence)<=1.5):
                diverg[j]=''
                for ind in indice:
                    diverg[j]=diverg[j]+subset.columns.values[ind]+'/'
            else:
                diverg[j]='_'
            diverg[j]=diverg[j][:-1]
    
    
    variant_id=list(np.arange(0,setargs.K))
    
    pd.DataFrame(diverg).to_csv(dir1+'diverg_id.csv',index=False)  
    true_theta['non-zero']=true_theta.astype(bool).sum(axis=1)
    true_theta['non-zero_prior']=true_theta.iloc[:,0:5].astype(bool).sum(axis=1)
    true_theta['non-zero_deno']=true_theta.iloc[:,5:7].astype(bool).sum(axis=1)
    
    if(setargs.K==2):
        mix=2;
        times=1;
        start=0;
    if(setargs.K==3):
        mix=3;
        times=1;
        start=0;
    if(setargs.K==7):
        mix=2;
        times=2;  
        start=0
    for i in range(0,times):    
        
        subset= true_theta[true_theta['non-zero']==mix]
        est_theta=pd.DataFrame(est_theta)
        subset_3=est_theta[est_theta.index.isin(subset.index)]
        
        subset_ref_ref=true_theta[(true_theta['non-zero_prior']==2) & (true_theta['non-zero_deno']==0) ]
        subset_4=est_theta[est_theta.index.isin(subset_ref_ref.index)]
        subset_ref_denovo=true_theta[(true_theta['non-zero_prior']==1) & (true_theta['non-zero_deno']==1)]
        subset_5=est_theta[est_theta.index.isin(subset_ref_denovo.index)]
        subset_ref_ref_denovo=true_theta[(true_theta['non-zero_prior']==2) & (true_theta['non-zero_deno']==1)]
        subset_6=est_theta[est_theta.index.isin(subset_ref_ref_denovo.index)]
        subset_ref_denovo_denovo=true_theta[(true_theta['non-zero_prior']==1) & (true_theta['non-zero_deno']==2)]
        subset_7=est_theta[est_theta.index.isin(subset_ref_denovo_denovo.index)]
        
        list1=[];
        list2=[];
        list3=[];
        list4=[];
        list5=[];
        list6=[];
        for k in range(0,setargs.K):
            list1.append(rela_error(subset.iloc[:,k],subset_3.iloc[:,int(diverg[k])]))
            list2.append(rela_error(true_theta.iloc[:,k],est_theta.iloc[:,int(diverg[k])]))
            list3.append(rela_error(subset_ref_ref.iloc[:,k],subset_4.iloc[:,int(diverg[k])]))
            list4.append(rela_error(subset_ref_denovo.iloc[:,k],subset_5.iloc[:,int(diverg[k])]))
            list5.append(rela_error(subset_ref_ref_denovo.iloc[:,k],subset_6.iloc[:,int(diverg[k])]))
            list6.append(rela_error(subset_ref_denovo_denovo.iloc[:,k],subset_7.iloc[:,int(diverg[k])]))
        
        rel_error_array[meth_count,start]=np.mean(list1);     #only mixed
        rel_error_array[meth_count,start+1]=np.mean(list1[0:prior]);
        rel_error_array[meth_count,start+2]=np.mean(list1[prior:setargs.K]);
        rel_error_array[meth_count,start+3]=np.mean(list2);   #all_samples including pure and mixed
        rel_error_array[meth_count,start+4]=np.mean(list3);   #all_samples including pure and mixed
        rel_error_array[meth_count,start+5]=np.mean(list4);   #all_samples including pure and mixed
        rel_error_array[meth_count,start+6]=np.mean(list5);   #all_samples including pure and mixed
        rel_error_array[meth_count,start+7]=np.mean(list6);   #all_samples including pure and mixed
        
        
        
        mix=3;
        start=8;
        
        with open(filename, 'wb') as f:
            pickle.dump(rel_error_array, f)
        f.close();
    
def runs(variants_dict,known,NMF=False,prior=False):
    var_initialize();
    
    print(setargs.K)
    
    if(NMF==False):
        print("running")
        os.system("../simulation1/Demixer 0 "+str(setargs.K))
    else:
        os.system("../simulation1/Demixer 1 "+str(setargs.K))
 
    filename = "Demixer_cpp/output/n_m_z.dat"
    #filename = "NIRT_Docs_int_onlySnpgap.dat"
    fileobj = open(filename, mode='r')
    setargs.n_m_z= np.fromfile(fileobj, dtype=np.uint32)
    fileobj.close
    
    
    
    filename = "Demixer_cpp/output/n_z_t.dat"
    #filename = "NIRT_Docs_int_onlySnpgap.dat"
    fileobj = open(filename, mode='r')
    setargs.n_z_t= np.fromfile(fileobj, dtype=np.uint32)
    fileobj.close

    print(setargs.n_m_z)
    setargs.n_m_z=setargs.n_m_z.reshape((len(setargs.Docs),setargs.K))
    setargs.n_z_t=setargs.n_z_t.reshape((setargs.K,setargs.V))
    
    est_theta=setargs.n_m_z/setargs.n_m_z.sum(axis=1)[:,None]
    
    
    if(NMF==False and prior==False):
            pd.DataFrame(est_theta).to_csv(dir1+'strain_'+str(no_strain)+'_SNP_LDA'+'.csv',index=False)
            pd.DataFrame(setargs.n_z_t).to_csv(dir1+'strain_'+str(no_strain)+'_SNP_LDA_NZT'+'.csv',index=False)
    elif(NMF==True and prior==False):
            pd.DataFrame(est_theta).to_csv(dir1+'strain_'+str(no_strain)+'_NMF_SNP_LDA'+'.csv',index=False)
            pd.DataFrame(setargs.n_z_t).to_csv(dir1+'strain_'+str(no_strain)+'_NMF_SNP_LDA_NZT'+'.csv',index=False)
    elif(NMF==False and prior==True):
            pd.DataFrame(est_theta).to_csv(dir1+'strain_'+str(no_strain)+'_Prior_SNP_LDA'+str(setargs.alpha)+'_'+str(setargs.beta)+'.csv',index=False)
            pd.DataFrame(setargs.n_z_t).to_csv(dir1+'strain_'+str(no_strain)+'_Prior_SNP_LDA_NZT'+str(setargs.alpha)+'_'+str(setargs.beta)+'.csv',index=False)
    elif(NMF==True and prior==True):
            pd.DataFrame(est_theta).to_csv(dir1+'strain_'+str(no_strain)+'_NMF_prior_SNP_LDA'+str(setargs.alpha)+'_'+str(setargs.beta)+'.csv',index=False)
            pd.DataFrame(setargs.n_z_t).to_csv(dir1+'strain_'+str(no_strain)+'_NMF_prior_SNP_LDA_NZT'+str(setargs.alpha)+'_'+str(setargs.beta)+'.csv',index=False)
    mutations(setargs.n_z_t,est_theta)            
    
def topic_kldiverg(no_strain,Docs_col_name,n_z_t,est_theta,filename):
    #print(Docs_col_name)
    K=setargs.K;
   
    true_theta=pd.read_csv('artsim/'+str(no_strain)+'strain_true_theta.csv')
    
    for i in range(0,setargs.K):
        true_theta[str(i)]=true_theta[(str(i))]/100
    
    new_topic=np.zeros(setargs.K);
    
    diverg=np.zeros(K)
    divergence=np.zeros(K);
    n_z_t=n_z_t+0.01
    for topic in range(0,K):
        divergence=np.zeros(K);
        divergence=divergence+diverg
        for i in range(0,int(len(Docs_col_name)/4)):
            if(Docs_col_name[i*4+1] in mut_dict):
           
                    id=mut_index[Docs_col_name[i*4+1]]
                    n_z_t_sub=n_z_t[:,i*4:(i+1)*4]
                    
                    n_z_sub=numpy.transpose(numpy.sum(n_z_t_sub,axis=1))
                    sub=n_z_t_sub / n_z_sub[:, numpy.newaxis]
                    sub=numpy.nan_to_num(sub, nan=0)
                    print(phi[id][topic][0],np.around(t,4))
                    divergence=divergence+np.array([entropy(phi[id][topic][0],np.around(t,4)) for t in sub])
        print(np.nanmin(divergence,np.where(divergence==np.nanmin(divergence))[0][0]))
        
        if(np.nanmin(divergence)<=0.5):
            indice=np.where(divergence==np.nanmin(divergence))[0][0]
            diverg[indice]=np.nan
            new_topic[topic]=indice
    indice=np.where(diverg==np.nanmin(diverg))[0][0]
    #mae_val[meth_count,topic+1]=rela_error(true_theta[:,topic+1], est_theta[:,indice])
    #new_topic[topic+1]=indice
    print(new_topic)
    true_theta['non-zero']=true_theta.astype(bool).sum(axis=1)
    
    
    if(setargs.K==2):
        mix=2;
        times=1;
        start=0;
    if(setargs.K==3):
        mix=3;
        times=1;
        start=0;
    if(setargs.K==7):
        mix=2;
        times=2;  
        start=0
    for i in range(0,times):    
        
        subset= true_theta[true_theta['non-zero']==mix]
        est_theta=pd.DataFrame(est_theta)
        subset_3=est_theta[est_theta.index.isin(subset.index)]
        
        
        list1=[];
        list2=[];
        for k in range(0,setargs.K):
            list1.append(rela_error(subset.iloc[:,k],subset_3.iloc[:,int(new_topic[k])]))
            list2.append(rela_error(true_theta.iloc[:,k],est_theta.iloc[:,int(new_topic[k])]))
        rel_error_array[meth_count,start]=np.mean(list1);     #only mixed
        rel_error_array[meth_count,start+1]=np.mean(list1[0:prior]);
        rel_error_array[meth_count,start+2]=np.mean(list1[prior:setargs.K]);
        rel_error_array[meth_count,start+3]=np.mean(list2);   #all_samples including pure and mixed
        rel_error_array[meth_count,start+4]=np.mean(list3);
        rel_error_array[meth_count,start+5]=np.mean(list4);
        rel_error_array[meth_count,start+6]=np.mean(list5);
        rel_error_array[meth_count,start+7]=np.mean(list6);
        mix=3;
        start=4;
        
        with open(filename, 'wb') as f:
            pickle.dump(rel_error_array, f)
        f.close();
    
def add_dict(s1,key,s2,s3,s4,s5,s6,s7):
    global variants_subset;
    v='mut'+s1
    v1='mut'+s2
    mutations=globals()[v]
    mut=globals()[v1]
    l1=[]
    for j in range(0,int(len(mutations))):
        c1=mutations[j];
        
        if(c1 in setargs.cols_name and c1 not in mut and c1 not in globals()['mut'+s3] and c1 not in globals()['mut'+s4] and c1 not in globals()['mut'+s5] and c1 not in globals()['mut'+s6] and c1 not in globals()['mut'+s7]):
                c=setargs.cols_name.get(c1)
                l1.append(c)
                if(c not in variants_dict):
                    variants_dict.setdefault(c, [])
                    variants_dict[int(c)].append(key);
    additional=pd.DataFrame()
    additional[key]=l1;
    variants_subset = pd.concat([variants_subset, additional], axis=1) 
    
def output_to_file(variants_dict):
    f = open('Demixer_cpp/input/temp1.txt', "w")
    for j, key in enumerate(variants_dict):
            
            l=variants_dict.get(key)
            if(type(l)==type([])):
                for i in l:
                    s1=""
                    #if(len(l)>0):
                    s1 = " ".join(str(i))
                    f.write(str(int(key))+ ' 1 ' +s1+'\n')
            else:
               s1=""
               s1 = " ".join(str(l))
               f.write(str(int(key))+ ' 1 ' +s1+'\n')
    f.close()                   
def add_dict1(s1,key,s2,s3,s4,s5):
    global variants_subset;
    v='mut'+s1
    v1='mut'+s2
    mutations=globals()[v]
    mut=globals()[v1]
    l1=[]
    for j in range(0,int(len(mutations))):
        c1=mutations[j];
        
        if(c1 in setargs.cols_name and c1 not in mut and c1 not in globals()['mut'+s3] and c1 not in globals()['mut'+s4] and c1 not in globals()['mut'+s5]):
                c=setargs.cols_name.get(c1)
                l1.append(c)
                if(c not in variants_dict):
                    variants_dict.setdefault(c, [])
                    variants_dict[int(c)].append(key);
    additional=pd.DataFrame()
    additional[key]=l1;
    variants_subset = pd.concat([variants_subset, additional], axis=1) 
    
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser

import sys; 
import os;

fasta_sequences = SeqIO.parse(open('../../scripts/db/tuberculosis.fna'),'fasta')
ref=str();
count=0;
total_len=0;
ref=[];
setargs.alpha=1;
setargs.beta=0.01;
no_strain=int(sys.argv[1]);


with open("../../scripts/db/tuberculosis.fna") as in_handle:
     for title, seq in SimpleFastaParser(in_handle):
            ref=seq
ch='A'

mut_dict=set();
mut_list=[];
phi=[];

if(no_strain==9):
        setargs.K=7
else:    
        setargs.K=no_strain;

prior=int(sys.argv[3]);

if(setargs.K==prior):
    dir1='artsim_errorbar'+sys.argv[2]+'/'
else:
    dir1='artsim_errorbar'+sys.argv[2]+'/denovo/'
    
isExist = os.path.exists(dir1)
if not isExist:

   # Create a new directory because it does not exist
    os.makedirs(dir1)
    print("The new directory is created!")

file_path = dir1+'KLoutput_strain'+str(no_strain)+'_'+sys.argv[4]+'timecheck.txt'

sys.stdout = open(file_path, "w")
        
no_mutations=np.zeros(setargs.K);

for j in range(0,setargs.K):           
    with open('artsim/'+str(no_strain)+'strain/'+str(no_strain)+'_strain'+ch+'_100.fasta') as in_handle: 
        for title, seq in SimpleFastaParser(in_handle):
            globals()[ch]=seq
    v='mut'+ch
    globals()[v]=[]
    
    for i in range(0, len(ref)):
        phi_q=[];
        if(ref[i]!=globals()[ch][i]):
            m=ref[i]+str(i+1)+globals()[ch][i]
            globals()[v].append(m)
            mut_list.append(m)
            mut_dict.add(m)
            for k in range(0,setargs.K):  
                if(j==k):
                    prob=[[0,1,0,0]];
                else:
                    prob=[[1,0,0,0]]
                #print(prob)
                phi_q.append(prob);
            phi.append(phi_q)
    #print(phi)
    ch=chr(ord(ch) + 1)
    
    

mut_index=dict(zip(mut_list,range(len(mut_list))))    

#setargs.K=int(sys.argv[2])
for i in range(0,1):
    variants_subset=pd.DataFrame()
    meth_count=0;

    filename='artsim/'+str(no_strain)+'strain.vcf'
    genotype=allel.vcf_to_dataframe(filename,fields='*'); 
    #genotype=cal_genotype(filename);
    print('reading done');
    #print(genotype)
    genotype=genotype[genotype['is_snp']==True]
    genotype=genotype.reset_index();
    row,col=genotype.shape 
    print(row,col);
    setargs.callset = allel.read_vcf(filename,fields='*')
    read_countdata=pd.DataFrame()
    #print(callset)
    read_countdata_proportion=pd.DataFrame()
    read_countdata_entropy=pd.DataFrame()
    sample_name=list()
    length=len(setargs.callset['samples'])
    #length=150;
    doc_n=length
    for j in range(0,length):     
            sample_name.append(setargs.callset['samples'][j]);
    read_countdata['Sample_id']=sample_name;
    temp=1;    

    
    print("Number of samples: "+ str(length))
    setargs.Docs=[[] for x in range(length)]
    setargs.Docs_1=[[] for x in range(length)]
    setargs.Docs_col_name=[];

    id_num=0
    for variant in range(0,row): 
        id_num=read_data(variant,id_num,genotype,length,read_countdata)  
    setargs.V=len(setargs.Docs_col_name);   #vocabulary
    setargs.cols_name=dict(zip(setargs.Docs_col_name,range(len(setargs.Docs_col_name))))
    setargs.idf= numpy.zeros(setargs.V)
    setargs.idf=setargs.idf.astype("int16")
   
    filepointer = open("Demixer_cpp/input/idf.dat", mode='wb')
    setargs.idf.tofile(filepointer)
    filepointer.close    
    
    with open('Demixer_cpp/input/config.txt', 'w') as file:
        file.write(str(len(setargs.Docs)) + ' '+str(setargs.K)+ ' '+str(setargs.V)+ ' '+str(setargs.V));
    
    Docs=np.array(setargs.Docs_1).astype("int16")

    setargs.Docs_1=np.array(setargs.Docs_1).astype("int16")

    filename = "Demixer_cpp/input/Docs_int_onlySnpgap.dat"
    fileobj = open(filename, mode='wb')
    Docs.tofile(fileobj)
    fileobj.close
    
    initial_dict=[];

    #parallel_initialize(initial_dict,prior,True);
    output_to_file({})
    parallel_construct(False,'Demixer_cpp/input/',0,3,pd.DataFrame())
    
    variants_dict={}
    ch='A';
    
    
    if(prior<5):
        id=0;
        for k in range(0,prior):
            l1=globals()['mut'+ch]
            l2=[]
            for j in range(0,len(l1)):
                c1=setargs.cols_name.get(l1[j])
                
                if(isinstance(c1, int)):
                    l2.append(c1)
                    if(c1 not in variants_dict):
                        variants_dict.setdefault(c1, [])
                        variants_dict[int(c1)].append(id);
            ch=chr(ord(ch)+1)
            id=id+1;
            additional=pd.DataFrame()
            additional[k]=l2;
            variants_subset = pd.concat([variants_subset, additional], axis=1) 

    elif(prior==5):

        add_dict('A',0,'B','C','D','E','F','G')
        add_dict('B',1,'A','C','D','E','F','G')
        add_dict('C',2,'B','A','D','E','F','G')
        add_dict('D',3,'B','C','A','E','F','G')
        add_dict('E',4,'B','C','D','F','G','A')

    else:
        add_dict('A',0,'B','C','D','E','F','G')
        add_dict('B',1,'A','C','D','E','F','G')
        add_dict('C',2,'B','A','D','E','F','G')
        add_dict('D',3,'B','C','A','E','F','G')
        add_dict('E',4,'B','C','D','F','G','A')
        add_dict('F',5,'B','C','D','E','G','A')
        add_dict('G',6,'B','C','D','E','F','A')
     
    '''
    df = pd.read_csv('artsim_errorbar/strain_'+str(no_strain)+'_Prior_SNP_LDA_NZT1_0.01.csv', sep=',', header=None)
    setargs.n_z_t=df.values[1:,:]
    df = pd.read_csv('artsim_errorbar/strain_'+str(no_strain)+'_Prior_SNP_LDA1_0.01.csv', sep=',', header=None)
    est_theta=df.values[1:,:]
    topic_kldiverg(setargs.K,setargs.Docs_col_name,setargs.n_z_t,est_theta,'artsim_errorbar/rel_error'+str(no_strain)+'.pkl')
    
    print(mutG)
   
    
    '''
    output_to_file(variants_dict)

    weight=parallel_construct(False,'Demixer_cpp/input/',0,3,pd.DataFrame(variants_subset))[1]
    print(weight)
   
    setargs.idf[variants_subset.stack().values.astype(int)]=weight    
    write_datfile("Demixer_cpp/input/idf.dat",setargs.idf)
    pd.DataFrame(variants_subset).to_csv('Demixer_cpp/input/subset.csv',index=False) 
    
    alpha=np.full(setargs.K,1,dtype='<f4')
    beta=np.full(setargs.K,0.01,dtype='<f4')
    write_datfile("Demixer_cpp/input/alpha.dat",alpha)
    write_datfile("Demixer_cpp/input/beta.dat",beta)
    
 
    
    j=int(sys.argv[4])
            
    print('running_prior_SNP_LDA',setargs.alpha,setargs.beta)



    if(no_strain==7 or no_strain==9):
        rel_error_array=np.zeros((1,16));
    else:
        rel_error_array=np.zeros((1,4));
    st = time.time()
    runs(variants_dict,prior,False,True)
    iter_time = time.time() - st
    print(iter_time)

    with open(dir1+'rel_error'+str(no_strain)+'_'+str(j)+'.pkl', 'wb') as f:
        pickle.dump(rel_error_array, f)
    f.close();

