import re,sys
import shutil
import os
import math
from operator import add 
import operator
import threading
import time
import numpy
from pandas import ExcelWriter
from pandas import ExcelFile
import pickle
import pdb;
import random
from sim_init import*
from nmf import*
from scipy.stats import entropy
import pandas as pd;
import numpy as np;
from multiprocessing.pool import ThreadPool
import multiprocessing as mp
import multiprocessing
import matplotlib.pyplot as plt
from plotnine import *
from scipy import stats
from sklearn.metrics import mean_absolute_error as mae
import pdb
import gc
current_directory = os.getcwd()
two_levels_up = os.path.abspath(os.path.join(current_directory, "../../scripts"))
sys.path.append(two_levels_up)

from setargs import *
from init import *
weight=0
setargs.strain=''
no_run=10
rel_error=np.zeros((3,6,no_run))
current_directory = os.getcwd()

coun=1
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
    

def simulate_data_3(K,SNP_NUM,READ_NUM,DOC_NUM,variants_dict):

    TOPIC_N = K;
    VOCABULARY_SIZE = SNP_NUM*4
    V=SNP_NUM*4;
    Docs=[];
    a_m_n = []
    a_m_z = numpy.zeros((DOC_NUM, K)) 
    a_z_t = numpy.zeros((K, V))
    a_z = numpy.zeros(K)
    a_m_t = numpy.zeros((DOC_NUM, V,K))
    a_r=numpy.zeros((DOC_NUM, int(V/4)+1))
    col_name=[]

    for i in range(SNP_NUM):
        for j in range(4):
            col_name.append(str(i)+'_'+str(j));
    #print(col_name);
    theta=[];
    start=100
    seq_error=[];

    for i in range(10):
        #print(start)
        for j in range(0,10):
            a1=1;
            a2=1;

            while(1-(a1+a2)<0):
                a= random.randint(0,20);
                a3=np.round((rng.normal(a))/100,2)
                a1=np.round((rng.normal(0.6,0.05)),2)
                #a2=np.round((rng.normal(start*40/100,10))/100,2)
                #a1=np.round((rng.normal(start,10))/100,2)
                #a2=np.round((rng.normal(start,10))/100,2)
                a2=np.round((1-(a1+a3)),3)
            #topic=[[start*60/10000,start*40/10000,(100-start)/100]]
           
            topic=[[a1,a2,a3]]
            #print(topic)
            theta.append(topic);
        start=start-5
    #print('Completed')
    #print(theta); 

    alpha_options=[(5,3,1),(5,3,1.5)]
    #rint(alpha_options[0])

    beta = numpy.full((SNP_NUM,3),0.01)
    #alpha = numpy.full(TOPIC_N,0.9);
    true_theta=numpy.zeros((DOC_NUM,TOPIC_N))
    #rint(beta,alpha)
 
    Dataset=pd.DataFrame(columns=col_name)
    phi = []
    ## generate multinomial distribution over words for each topic
    for i in range(SNP_NUM):
        phi_q=[];
        flag=0;
        if i in set(variants.iloc[:,0]):
          
            topic=[[0,1,0,0]];
            phi_q.append(topic);
            topic=[[1,0,0,0]];
            phi_q.append(topic);

            topic=[[1,0,0,0]];
            phi_q.append(topic);


        elif i in set(variants.iloc[:,1]):
            topic=[[1,0,0,0]];
            phi_q.append(topic);
            topic=[[0,1,0,0]];
            phi_q.append(topic);

            topic=[[1,0,0,0]];
            phi_q.append(topic);

        elif i in set(variants.iloc[:,2]):

            topic=[[1,0,0,0]];
            phi_q.append(topic);
            topic=[[1,0,0,0]];
            phi_q.append(topic);

            topic=[[0,1,0,0]];
            phi_q.append(topic);


        phi.append(phi_q);
    #print(phi)
    ## generate words for each document
    words=np.zeros(V);
    start=100;
    for i in range(DOC_NUM): 
        z_n = []
        t1=0;
        terms=[];
        buffer = {}
        z_buffer = {} ## keep track the true z
        ## first sample theta
        #print(i)
        #theta = rng.mtrand.dirichlet(alpha_options[rng.randint(0, 2)],size = 1)
        #print(theta[i])
        row_list=np.zeros(V);
        for j in range(SNP_NUM):
                for jj in range(READ_NUM):
                
                    z = rng.multinomial(1,theta[i][0],size = 1)
                    #print('z',j,z);
                    z_assignment = 0
                    for k in range(TOPIC_N):
                        if z[0][k] == 1:
                            break
                        z_assignment += 1
                    if not z_assignment in z_buffer:
                        z_buffer[z_assignment] = 0
                    z_buffer[z_assignment] = z_buffer[z_assignment] + 1
                     
                    w = rng.multinomial(1,phi[j][z_assignment][0],size = 1)
                    #print(i,theta[i][0],z_assignment,'i=',i,'j=',j,phi[j],"w=",w)
                    #if(j==100):
                    #print(z_assignment);
                    #print('w',w
                    w_assignment = 0
                    for k in range(4):
                        if w[0][k] == 1:
                            break
                        w_assignment += 1
                    if not w_assignment in buffer:
                        buffer[w_assignment] = 0
                    buffer[w_assignment] = buffer[w_assignment] + 1

                    z_n.append(z)
                    r=j*4+w_assignment
                    a_m_z[i, z_assignment] += 1
                    a_m_t[i,r,z_assignment]+= 1
                    a_z_t[z_assignment, r] += 1
                    a_z[z_assignment] += 1
                    t1+=1;
                    rr=int(r/4)
                    a_r[i,rr+1]+=1
                    #n_r_max[m,rr+1]=max(n_r_max[m,rr],na)
                    row_list[j*4+w_assignment]+=1;
                    terms.append(j*4+w_assignment);
                    setargs.Docs_1[i][j*4+w_assignment]+=1
        #print(row_list)
        Dataset.loc[len(Dataset)] = row_list
        a_m_n.append(numpy.array(z_n))
        Docs.append(terms)
        #print(Dataset)
      
        true_theta[i]=theta[i][0]
       
 
    with open('synthetic'+str(no_strain)+'strain/phi.pkl', 'wb') as f:
                pickle.dump(phi, f)
    f.close();
    
    return Docs, true_theta, Dataset

def topic_kldiverg(K,Docs_col_name,n_z_t,est_theta):
    global coun;
    #print(Docs_col_name)
    mut_list=[];
    for j in range(0,K):
        mut_list.append([]);
    unpickleFile = open('synthetic'+str(no_strain)+'strain/phi.pkl', 'rb')
    phi = pickle.load(unpickleFile)
    
    diverg=np.zeros(K)
    divergence=np.zeros(K);
    
    n_z_t=n_z_t+0.01
    n_z_t_prop=np.zeros((K,len(Docs_col_name)*4))
        
    for i in range(0,int(len(Docs_col_name)/4)):
        n_z_t_sub=n_z_t[:,i*4:(i+1)*4]
        n_z_sub=numpy.transpose(numpy.sum(n_z_t_sub,axis=1))
        sub=n_z_t_sub / n_z_sub[:, numpy.newaxis]
        sub=numpy.nan_to_num(sub, nan=0)      
        n_z_t_prop[:,i*4:(i+1)*4] = sub
    
    for topic in range(0,K):
        divergence=np.zeros(no_strain);
        for i in range(0,int(len(Docs_col_name)/4)):
            sub=n_z_t_prop[:,i*4:(i+1)*4]
            sub=sub+0.01
            for j in range(0,no_strain):
                    phi_ext=[0,0,0,0]
                    for jj in range(0,4):
                        
                        phi_ext[jj]=phi[i][j][0][jj]+0.01
      
                    #print(np.around(sub[topic,:],4),sub[topic,:],phi_ext)
                    divergence[j]=divergence[j]+np.array([entropy(phi[i][j][0],np.around(sub[topic,:],4))])
        
        
        divergence=divergence/(len(Docs_col_name)/4)
   
        if(np.nanmin(divergence)<=1.5):     
            indice=np.where(divergence==np.nanmin(divergence))[0][0]
            #mae_val[meth_count,topic]=rela_error(true_theta[:,topic], est_theta[:,indice])
            diverg[topic]=indice
            print(divergence)
            #diverg[indice]=np.nan
        if(np.nanmin(divergence)>1.5):
                print(divergence)
                print("greater than 1.5")
                
    print(diverg)
    
    c_theta=pd.DataFrame()
    for topic in range(0,no_strain):
        indice=np.where(diverg==topic)
        prop=np.zeros(len(true_theta[:,topic]))
        if(len(indice[0]>0)):
            for ind in indice[0]:
                prop=prop + est_theta[:,ind]
        mae_val[meth_count,topic]=rela_error(true_theta[:,topic], prop)
        c_theta[topic]=prop
        
    print("meth",meth_count)
    
    pd.DataFrame(est_theta).to_csv(r'synthetic'+str(no_strain)+'strain/demixer_'+str(meth_count)+'_'+str(coun)+'.csv',index=False)
    coun+=1
    pd.DataFrame(c_theta).to_csv(r'synthetic'+str(no_strain)+'strain/demixer_'+str(meth_count)+'_'+str(coun)+'.csv',index=False)
    
    coun+=1
    print(topic+1,indice)
    
def generate_SNPs(no_strains,no_snps):
    variants=pd.DataFrame()
    variants_dict=dict();
    l=list(np.arange(no_strains*no_snps))
    
    random.shuffle(l);
    
    #print(np.shuffle(l))
    for i in range(0,no_strains):
        l1=l[i*no_snps:(i+1)*no_snps];
        for j in range(0,len(l1)):
            c1=l1[j]*4 +1;
            if(c1 not in variants_dict):
                variants_dict[c1]=i;
        variants['Type_'+str(i)]=l1;
    
    return variants, variants_dict;
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
    
    
def simulate_data_4(K,SNP_NUM,READ_NUM,DOC_NUM,variants_dict):

    TOPIC_N = K;
    VOCABULARY_SIZE = SNP_NUM*4
    V=SNP_NUM*4;
    Docs=[];
    a_m_n = []
    a_m_z = numpy.zeros((DOC_NUM, K)) 
    a_z_t = numpy.zeros((K, V))
    a_z = numpy.zeros(K)
    a_m_t = numpy.zeros((DOC_NUM, V,K))
    a_r=numpy.zeros((DOC_NUM, int(V/4)+1))
    col_name=[]

    for i in range(SNP_NUM):
        for j in range(4):
            col_name.append(str(i)+'_'+str(j));
    theta=[];
   
    for i in range(10):
        for j in range(0,10):
            a1=1;
            a2=1;

            while(1-(a1+a2)<0):
                a= random.randint(0,20);
                a3=np.round((rng.normal(a))/100,2)
                a1=np.round((rng.normal(0.6,0.05)),2)
                a2=np.round((1-(a1+a3)),3)
            topic=[[a1,a2,a3,0]]
    
            theta.append(topic);


    #print(theta)
    for i in range(10):
        a=random.randint(0,99)
        a3=abs(np.round(rng.normal(0.1,0.05),2))
        theta[a][0][3]=a3;
        theta[a][0][0]= np.round(theta[a][0][0] - np.round(a3/2,4),2)
        theta[a][0][1]= np.round(theta[a][0][1] - np.round(a3/2,4),2)

    true_theta=numpy.zeros((DOC_NUM,TOPIC_N))

  
    Dataset=pd.DataFrame(columns=col_name)
    phi = []
    ## generate multinomial distribution over words for each topic
    for i in range(SNP_NUM):
        phi_q=[];
        flag=0;
        
        for j in range(0,4):
            if i in set(variants.iloc[:,j]):
                for k in range(0,K):
                    if(j==k):
                        topic=[[0,1,0,0]];
                    else:
                        topic=[[1,0,0,0]];
                    phi_q.append(topic);
        #print(phi_q)
        phi.append(phi_q);
            
    ## generate words for each document
    words=np.zeros(V);
    
    for i in range(DOC_NUM): 
        z_n = []
        t1=0;
        terms=[];
        buffer = {}
        z_buffer = {} ## keep track the true z
     
        row_list=np.zeros(V);
        for j in range(SNP_NUM):
                for jj in range(READ_NUM):
            ## first sample z
                    z = rng.multinomial(1,theta[i][0],size = 1)
                    #print('z',j,z);
                    z_assignment = 0
                    for k in range(TOPIC_N):
                        if z[0][k] == 1:
                            break
                        z_assignment += 1
                    if not z_assignment in z_buffer:
                        z_buffer[z_assignment] = 0
                    z_buffer[z_assignment] = z_buffer[z_assignment] + 1
            ## sample a word from topic z            
                    #print(phi[j][z_assignment])
                    w = rng.multinomial(1,phi[j][z_assignment][0],size = 1)
                    
                    w_assignment = 0
                    for k in range(4):
                        if w[0][k] == 1:
                            break
                        w_assignment += 1
                    if not w_assignment in buffer:
                        buffer[w_assignment] = 0
                    buffer[w_assignment] = buffer[w_assignment] + 1

                    z_n.append(z)
                    r=j*4+w_assignment
                    a_m_z[i, z_assignment] += 1
                    a_m_t[i,r,z_assignment]+= 1
                    a_z_t[z_assignment, r] += 1
                    a_z[z_assignment] += 1
                    t1+=1;
                    rr=int(r/4)
                    a_r[i,rr+1]+=1
                    #n_r_max[m,rr+1]=max(n_r_max[m,rr],na)
                    row_list[j*4+w_assignment]+=1;
                    terms.append(j*4+w_assignment);
                    setargs.Docs_1[i][j*4+w_assignment]+=1
        #print(row_list)
        #print(a_m_z)
        
        Dataset.loc[len(Dataset)] = row_list
        a_m_n.append(numpy.array(z_n))
        Docs.append(terms)
        #print(Dataset)
        true_theta[i]=theta[i][0]

    with open('synthetic'+str(no_strain)+'strain/phi.pkl', 'wb') as f:
                pickle.dump(phi, f)
    f.close();
    
    return Docs, true_theta, Dataset




def var_initialize():
    
    setargs.z_m_n=[];
    #setargs.strain=sys.argv[3];    #output_filename
    setargs.n_m_z = numpy.zeros((len(setargs.Docs), setargs.K),dtype='uint32')     # word count of each document and topic
    setargs.n_z_t = numpy.zeros((setargs.K, setargs.V),dtype='uint32') # word count of each topic and vocabulary
    setargs.n_z = numpy.zeros(setargs.K)        # word count of each topic
    setargs.n_m_t = numpy.zeros((len(setargs.Docs), setargs.V,setargs.K),dtype='int16')
    setargs.n_r=numpy.zeros((len(setargs.Docs), int(setargs.V/4)+1))
    setargs.n_r_max=numpy.full((len(setargs.Docs), int(setargs.V/4)+1),-1)
    setargs.doc_term= numpy.zeros((len(setargs.Docs), setargs.V))
    setargs.idf= numpy.zeros(setargs.V)

def runs(variants_dict,known,NMF=False,prior=False):
    var_initialize();
    
    
    ss_run=numpy.random.SeedSequence()
    #ss = rng.SeedSequence(entropy=289770586330054155256335069517751853071)
    print(ss_run.entropy)
    child_states_run = ss_run.spawn(300)
    if(NMF==False):
     
        #parallel_initialize(variants_dict,known,True);
        os.system("./Demixer 0 "+str(known))
    else:
        NMF_parallel_initialize(variants_dict,H,W);
        os.system("./Demixer 1 "+str(known))
 
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
    
    setargs.n_m_z=setargs.n_m_z.reshape((len(setargs.Docs),setargs.K))
    setargs.n_z_t=setargs.n_z_t.reshape((setargs.K,setargs.V))
    
    
    est_theta=setargs.n_m_z/setargs.n_m_z.sum(axis=1)[:,None]
    
    
    if(NMF==False and prior==False):
            pd.DataFrame(est_theta).to_csv(r'synthetic'+str(no_strain)+'strain/SNP_LDA'+str(run)+'.csv',index=False)
            pd.DataFrame(setargs.n_z_t).to_csv(r'synthetic'+str(no_strain)+'strain/SNP_LDA_NZT'+str(run)+'.csv',index=False)
    elif(NMF==True and prior==False):
            pd.DataFrame(est_theta).to_csv(r'synthetic'+str(no_strain)+'strain/NMF_SNP_LDA'+str(run)+'.csv',index=False)
            pd.DataFrame(setargs.n_z_t).to_csv(r'synthetic'+str(no_strain)+'strain/NMF_SNP_LDA_NZT'+str(run)+'.csv',index=False)
    elif(NMF==False and prior==True):
            pd.DataFrame(est_theta).to_csv(r'synthetic'+str(no_strain)+'strain/Prior_SNP_LDA'+str(run)+str(setargs.alpha)+'_'+str(setargs.beta)+'.csv',index=False)
            pd.DataFrame(setargs.n_z_t).to_csv(r'synthetic'+str(no_strain)+'strain/Prior_SNP_LDA_NZT'+str(run)+str(setargs.alpha)+'_'+str(setargs.beta)+'.csv',index=False)
    elif(NMF==True and prior==True):
            pd.DataFrame(est_theta).to_csv(r'synthetic'+str(no_strain)+'strain/NMF_prior_SNP_LDA'+str(run)+str(setargs.alpha)+'_'+str(setargs.beta)+'.csv',index=False)
            pd.DataFrame(setargs.n_z_t).to_csv(r'synthetic'+str(no_strain)+'strain/NMF_prior_SNP_LDA_NZT'+str(run)+str(setargs.alpha)+'_'+str(setargs.beta)+'.csv',index=False)
        
    topic_kldiverg(setargs.K,setargs.Docs_col_name,setargs.n_z_t,est_theta)

import sys; 
import os;
from numpy.random import SeedSequence, default_rng
setargs.alpha=1;
setargs.beta=0.01;
no_strain=int(sys.argv[1]);

path = os.path.join(current_directory, r'synthetic'+str(no_strain)+'strain/')


if(not os.path.exists(path)):
    os.mkdir(path)
#file_path = 'synthetic'+str(no_strain)+'strain/KLoutput_strain'+str(no_strain)+'.txt'
#sys.stdout = open(file_path, "w")

ss=numpy.random.SeedSequence()
#ss = rng.SeedSequence(entropy=289770586330054155256335069517751853071)
print(ss.entropy)
child_states = ss.spawn(no_run)

#try:
  
for i in range(0,no_run):
            #if(hk==0):
            setargs.K=10
            #else:
 

            rng=numpy.random.default_rng(child_states[i])
            meth_count=0;

            mae_val=np.zeros((6,setargs.K))
            run=i;
            
            
            for filename in os.listdir("Demixer_cpp/input/"):
                file_path = os.path.join("Demixer_cpp/input/", filename)
                if os.path.isfile(file_path):  # Check if it's a file (not a subdirectory)
                    os.remove(file_path) 

            variants,variants_dict=generate_SNPs(no_strain,50);
            #output_to_file(variants_dict)


            print('simulate data run=',i)
            setargs.Docs_1=[[0]*(no_strain*50*4) for x in range(0,100)]
            if(no_strain==3):
                setargs.Docs,true_theta, Dataset=simulate_data_3(no_strain,no_strain*50,50,100,variants_dict);
            else:
                setargs.Docs,true_theta, Dataset=simulate_data_4(no_strain,no_strain*50,50,100,variants_dict);

            pd.DataFrame(true_theta).to_csv(r'synthetic'+str(no_strain)+'strain/True_theta'+str(run)+'.csv',index=False)

            print(len(setargs.Docs[0]),true_theta.shape)


            Docs=np.array(setargs.Docs_1).astype("int16")

            setargs.Docs_1=np.array(setargs.Docs_1).astype("int16")


            write_datfile("Demixer_cpp/input/Docs_int_onlySnpgap.dat",Docs)
            #setargs.Docs=np.array(setargs.Docs)

            setargs.V=no_strain*50*4

            with open('Demixer_cpp/input/config.txt', 'w') as file:
                file.write(str(len(setargs.Docs)) + ' '+str(setargs.K)+ ' '+str(setargs.V)+ ' '+str(setargs.V)+ ' '+str(setargs.V));
            
            setargs.Docs_col_name=Dataset.columns[0:].values
            
            initial_dict=[];

            var_initialize();

            sim_parallel_initialize(initial_dict,setargs.K,True);

            
            output_to_file({})
            parallel_construct(False,'Demixer_cpp/input/',0,3,pd.DataFrame())
            
            #W,H=run_NMF(setargs.doc_term,no_strain,'synthetic'+str(no_strain)+'strain/',run)
            #est_theta=W/W.sum(axis=1)[:,None]


            #topic_kldiverg(setargs.K,setargs.Docs_col_name,H,est_theta)
            #no_mutations_l(setargs.K,setargs.Docs_col_name,H,variants)
            setargs.idf= numpy.zeros(setargs.V)
            setargs.idf=setargs.idf.astype("int16")

            ###################No reference high K ################################
            
            
            print('running_SNP_LDA')
            setargs.idf= numpy.zeros(setargs.V)
            setargs.idf=setargs.idf.astype("int16")
            write_datfile("Demixer_cpp/input/idf.dat",setargs.idf)
            
            alpha=np.full(setargs.K,1,dtype='<f4')
            beta=np.full(setargs.K,0.01,dtype='<f4')
            write_datfile("Demixer_cpp/input/alpha.dat",alpha)
            write_datfile("Demixer_cpp/input/beta.dat",beta)

            
            runs(initial_dict,0);

            
            gc.collect()
            
            
            meth_count+=1;
            ###################No reference Kset ################################
            
            output_to_file({})
            print('running_SNP_LDA')
            if(no_strain==3):
                setargs.K=3
            else:
                setargs.K=4
            
            with open('Demixer_cpp/input/config.txt', 'w') as file:
                file.write(str(len(setargs.Docs)) + ' '+str(setargs.K)+ ' '+str(setargs.V)+ ' '+str(setargs.V)+ ' '+str(setargs.V));


            runs(initial_dict,0);

            
            gc.collect()
            
            meth_count+=1;
            ###################reference + high K + wtset ################################
            #setargs.K=10
            variants_dict={}
            l1=variants['Type_0'].values
            l2=[]
            for j in range(0,len(l1)):
                c1=l1[j]*4 +1;
                l2.append(c1);
                if(c1 not in variants_dict):
                    variants_dict.setdefault(c1, [])
                    variants_dict[int(c1)].append(0);
                    
            output_to_file(variants_dict)
            sub_variants=pd.DataFrame()
            sub_variants['Type_0']=l2
            
            

            with open('Demixer_cpp/input/config.txt', 'w') as file:
                file.write(str(len(setargs.Docs)) + ' '+str(setargs.K)+ ' '+str(setargs.V)+ ' '+str(setargs.V)+ ' '+str(setargs.V));

            weight=parallel_construct(False,'Demixer_cpp/input/',0,3,sub_variants)[1]
            setargs.idf[l2]=weight
            write_datfile("Demixer_cpp/input/idf.dat",setargs.idf)


            print('running reference + highK + wtset',setargs.alpha,setargs.beta)

            runs(variants_dict,1,False,True)
            
            meth_count+=1;
            gc.collect()

            ###################reference + high K + wt not set ################################
            
            weight=0
            setargs.idf[l2]=weight
            write_datfile("Demixer_cpp/input/idf.dat",setargs.idf)

            print('running reference + highK + wt not set',setargs.alpha,setargs.beta)

            runs(variants_dict,1,False,True)
            
            meth_count+=1;
            gc.collect()
            ###################reference + K set + wt not set ################################
            if(no_strain==3):
                setargs.K=3
            else:
                setargs.K=4

            with open('Demixer_cpp/input/config.txt', 'w') as file:
                file.write(str(len(setargs.Docs)) + ' '+str(setargs.K)+ ' '+str(setargs.V)+ ' '+str(setargs.V)+ ' '+str(setargs.V));
            
            
            weight=0
            setargs.idf[l2]=weight
            write_datfile("Demixer_cpp/input/idf.dat",setargs.idf)

            print('running reference+ K+ wt not set',setargs.alpha,setargs.beta)

            
            runs(variants_dict,1,False,True)
            
            meth_count+=1;

            gc.collect()

            ###################reference + K set + wt set ################################
            weight=parallel_construct(False,'Demixer_cpp/input/',0,3,sub_variants)[1]
            setargs.idf[l2]=weight
            write_datfile("Demixer_cpp/input/idf.dat",setargs.idf)

            print('running reference + K set + wt set',setargs.alpha,setargs.beta)


            runs(variants_dict,1,False,True)
            gc.collect()

            for j in range(0,6):
                        rel_error[0,j,run]=np.mean(mae_val[j,:])
                        rel_error[1,j,run]=mae_val[j,0]
                        rel_error[2,j,run]=np.mean(mae_val[j,1:])


            with open('rel_error_'+str(no_strain)+'strains_ablation.pkl', 'wb') as f:
                pickle.dump(rel_error, f)
            f.close();
#except Exception as e:
#    print(f"An error occurred: {e}")
