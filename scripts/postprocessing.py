import numpy
import math
from scipy.stats import entropy
import numpy as np
import os
import pandas as pd
import multiprocessing
import pickle
import warnings
import matplotlib.pyplot as plt
import seaborn as sns
import math
import sys
import pdb
from scipy import stats
from setargs import *

warnings.filterwarnings("ignore")

def mutations(n_z_t,denovo_merge=True):
        
    subset=pd.read_csv(dir+'subset.csv')
    
    row,col=subset.shape;
    r,c=n_z_t.shape
    print(row,col)
    true_n_z_t = np.empty((col,c))
    true_n_z_t[:] = np.NaN
    l=len(subset.columns.values)
    d_count=0;
    for m in range(0,l):
            variant=list(subset.iloc[:,m])
            
            for i in range(0,len(variant)):
                
                    
                    if(not math.isnan(variant[i])):
                        variant[i]=int(variant[i])
                        ind=(variant[i]-snppos[variant[i]])%snptype[variant[i]]
                        #if(subset.columns[m]!='4' and subset.columns[m]!='4.9'):
                        loc=variant[i]-ind
                        #print(loc,snptype[loc],ind)
                        if(ind==0):
                            true_n_z_t[m,loc:loc+snptype[loc]]=np.zeros(snptype[loc])
                        elif(ind==1):
                            true_n_z_t[m,loc:loc+snptype[loc]]=np.zeros(snptype[loc])
                        elif(ind==2):
                            true_n_z_t[m,loc:loc+snptype[loc]]=np.zeros(snptype[loc])
                        elif(ind==3):
                            true_n_z_t[m,loc:loc+snptype[loc]]=np.zeros(snptype[loc])
                            
                        true_n_z_t[m,variant[i]]=1
    
    diverg=["" for x in range(0,K-len(exclude))]
    incre=0
    for j in range(K):
        if(j not in exclude):
            divergence=np.zeros(l);
            for m in range(0,l):
                variant=list(subset.iloc[:,m])
                count=0;
                for i in range(0,len(variant)):
                    if(not math.isnan(variant[i])):
                            #print(i,variant[i])
                            variant[i]=int(variant[i])
                            
                            ind=snptype[variant[i]]
                            loc=variant[i]-(variant[i]-snppos[variant[i]])%snptype[variant[i]]
                            
                            a=true_n_z_t[m][loc:loc+ind]
                            b=phi1[j][loc:loc+ind]
                            
                            divergence[m]=divergence[m]+entropy(a,b)
                            count+=1
                            if(j==2):
                                print(j,m,count,a,b)
                divergence[m]=np.around(divergence[m]/count,2)
                
            indice=np.where(divergence==np.nanmin(divergence))[0]
            print(j,indice,np.nanmin(divergence))
            
            #if(not math.isinf(np.nanmin(divergence))):
            if(np.nanmin(divergence)<=1.5):
                diverg[incre]=''
                #for ind in indice:
                diverg[incre]=subset.columns.values[indice[-1]]
            else:
                if(denovo_merge==True):
                    diverg[incre]='d_'
                else:
                    diverg[incre]='d_'+str(d_count)
                    d_count+=1
            #diverg[incre]=diverg[incre][:-1]
            incre+=1
    
    variant_id=list(np.arange(0,K))
    
    pd.DataFrame(diverg).to_csv(dir+'diverg_id.csv',index=False)  
    pd.DataFrame(theta,columns=diverg).to_csv(dir+'Prior.csv',index=False)
    
    return diverg

def mutations_fetch(n_z_t):
    
    for j in range(0,K):
        #print(j)
        mut=[];
        for i in range(0,V):
            if((i-snppos[i])%snptype[i]==0):
                #if(np.max(phi1[j,i:i+snptype[i]])!=0 and np.argmax(phi1[j,i:i+snptype[i]])!=0):
                mut.append(Docs_col_name[i+np.argmax(phi1[j,i:i+snptype[i]])])
        
        #print(len(mut))
        mut_l.append(mut)   

def kldiverg_variation(topic_num,K,V,index,exclude):
    #print(topic_num)
    topic=topic_num
    #print(K,V)
    start=0;
    
    for i in range(0,K):
        if(i not in exclude):
            ii=start
            if(index!=ii):
                    Vcount=0;
                    for j in range(0,V):
                        if((j-snppos[j])%snptype[j]==0):
                            temp=entropy(phi1[topic,j:j+snptype[j]], phi1[i,j:j+snptype[j]])
                            #if(not math.isnan(temp)):   
                            Vcount+=1
                            diverg_m[index,ii]=diverg_m[index,ii]+temp
                    diverg_m[index,ii]=(diverg_m[index,ii])/(Vcount)
            else:
                diverg_m[index,ii]=0
            start+=1


for i in range(1,2):
    #dir="../scripts/finaloutput/Cryptic_1st_run/"
    #dir="../scripts/finaloutput/Malawi_reset/"
    
    dir=sys.argv[1]
    print(dir)
    #os.chdir(dir)
    file_path = dir+"config.txt"
    filename = dir+"Docs_int_onlySnpgap.dat"
    fileobj = open(filename, mode='r')
    Docs_1 = np.fromfile(fileobj, dtype=np.int16)
    fileobj.close
    with open(dir+'Docs_col_name.pkl', 'rb') as f:
        Docs_col_name=pickle.load(f)
    f.close();
   
    with open(file_path, "r") as file:
        line = file.readline()
        integer_strings = line.split()
        integers = [int(i) for i in integer_strings]
    print(integers)
    m=integers[0]
    K=integers[1]
    V=integers[2]
    snp1=integers[3]
    snp2=integers[4]
   

    filename = dir+"n_m_z0.dat"
    Docs_1=Docs_1.reshape(m,V)
    fileobj = open(filename, mode='r')
    n_m_z = np.fromfile(fileobj, dtype=np.uint32)
    fileobj.close
    n_m_z=n_m_z.reshape((m,K))
    
    theta=n_m_z/n_m_z.sum(axis=1)[:,None]

    theta=np.around(theta,2)
    theta=theta/theta.sum(axis=1)[:,None]
    theta=np.around(theta,2)
    filename = dir+"n_z_t0.dat"
    fileobj = open(filename, mode='r')
    n_z_t = np.fromfile(fileobj, dtype=np.uint32)
    fileobj.close
    n_z_t=n_z_t.reshape((K,V))
    #n_z_t=n_z_t[0:15,:]



    snptype=np.zeros(V,dtype='int32');
    snppos=np.zeros(V,dtype='int32');


    snptype=multiprocessing.Array('i',V, lock=False)
    snptype = np.ctypeslib.as_array(snptype)

    snppos=multiprocessing.Array('i',V, lock=False)
    snppos = np.ctypeslib.as_array(snppos)

    print(n_z_t.shape)


    for j in range(0,V):
        if(j<integers[3]):
            snptype[j]=2;
            snppos[j]=0;
        elif(j<integers[4]):
            snptype[j]=3;
            snppos[j]=integers[3];
        else:
            snptype[j]=4;
            snppos[j]=integers[4];
    exclude=[]
    #phi1=np.zeros((K,V))
    for i in range(0,K):
        l1=list(np.where(theta[:,i]>0))
        if(len(l1[0])==0):
            print(i,l1)
            exclude.append(i)
    #exclude=[]
    theta=np.delete(theta, exclude, 1)
    print("exclude", exclude)
    pd.DataFrame(theta).to_csv(dir+'Prior.csv',index=False)

    l=len(exclude)
    diverg_m=multiprocessing.Array('d',(K-l)*(K-l), lock=False)
    diverg_m = np.ctypeslib.as_array(diverg_m)
    diverg_m = diverg_m.reshape(K-l,K-l)

    phi1=multiprocessing.Array('d',K*V, lock=False)
    phi1 = np.ctypeslib.as_array(phi1)
    phi1 = phi1.reshape(K,V)
    n_z_t=n_z_t+0.01
    for i in range(0,V):

            if((i-snppos[i])%snptype[i]==0):
                s_u_m=n_z_t[:,i:i+snptype[i]].sum(axis=1)[:,None]
                n_z_t_sub=n_z_t[:,i:i+snptype[i]]/s_u_m
                phi1[:,i:i+snptype[i]]=numpy.nan_to_num(n_z_t_sub, nan=0)    

    #phi=phi1+0.01
    '''
    process_pool = multiprocessing.Pool(processes=100)
    data=[]
    start1=0
    for k in range(0,K):
        if(k not in exclude):
            data.append([k,K,V,start1,exclude]);
            start1=start1+1;
    process_pool.starmap(kldiverg_variation, data)
    process_pool.close()

    pd.DataFrame(np.around(diverg_m,3)).to_csv(dir+'kldiverg.csv',index=False) 
    
    '''
    if(len(sys.argv)==3):
        diverg=mutations(n_z_t,sys.argv[2])
    else:
        diverg=mutations(n_z_t)
    

data=pd.DataFrame(Docs_1,columns=Docs_col_name)
cols_name=dict(zip(Docs_col_name,range(len(Docs_col_name))))    
'''
mut_l=[]
mutations_fetch(n_z_t)
variants=pd.read_csv('../scripts/db/tbprof_snpset.csv')        #for invitro
    
seeds=variants.stack().values 

mut_list=mut_l.copy();

data=[]
start1=0
for i in range(0,K):
    for j in range(0,K):
            if(i!=j):
                unwanted_num=list(set(mut_l[i]) & set(mut_l[j]))
                mut_list[i]=list(set(mut_list[i]) - set(unwanted_num))
              
'''
t=[]
long_ind=[]
sam_mean=[]
#Docs_col_name=[]
for i in range(0,m):
        #if(Demixer_result[0][i]=="Mix"):
        mode,l,mean=plot_prop_major(Docs_1[i],snp1,snp2,len(Docs_col_name))
        t.append(mode)
        long_ind.append(l)
        sam_mean.append(mean)
snp_plot_conf=pd.DataFrame()
snp_plot_conf['avg']=t
snp_plot_conf['long_bar']=long_ind
snp_plot_conf['mean']=sam_mean
snp_plot_conf.to_csv(dir+'snp_plot_conf.csv',index=False)

thresh_conf(dir)


