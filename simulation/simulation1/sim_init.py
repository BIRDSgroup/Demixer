import numpy
import numpy as np;
from multiprocessing import set_start_method
import multiprocessing
import pandas as pd
import random
import warnings
import math
import random
import os
import sys

current_directory = os.getcwd()
two_levels_up = os.path.abspath(os.path.join(current_directory, "../../scripts"))
sys.path.append(two_levels_up)

from setargs import *
import setargs

warnings.filterwarnings("ignore")
os.environ['OPENBLAS_NUM_THREADS'] = '1'

def child_initialize(_doc,variants_dict_,n_m_z_,n_m_t_,n_r_,doc_term_,known_,synthetic_,K_):
    global docs, _variants_dict, _n_m_z,_n_m_t,_n_r,_doc_term,known,synthetic,K
    docs= _doc
    _variants_dict= variants_dict_
    _n_m_z=n_m_z_
    _n_m_t=n_m_t_
    _n_r=n_r_
    _doc_term=doc_term_
    known=known_
    synthetic=synthetic_
    K=K_
    
def NMF_child_initialize(_doc,variants_dict_,n_m_z_,n_m_t_,n_r_,doc_term_,theta_est_,word_topic_,K_,V_):
    global docs, variants_dict, _n_m_z,_n_m_t,_n_r,_doc_term,_theta_est,word_topic,K,V
    docs= _doc
    variants_dict= variants_dict_
    _n_m_z=n_m_z_
    _n_m_t=n_m_t_
    _n_r=n_r_
    _doc_term=doc_term_
    K=K_
    V=V_
    _theta_est=theta_est_
    word_topic=word_topic_
    
def sim_parallel_initialize(variants_dict,known,synthetic=False):

        process_pool = multiprocessing.Pool(100,initializer = child_initialize, initargs = (setargs.Docs,variants_dict,setargs.n_m_z,setargs.n_m_t,setargs.n_r,setargs.doc_term,known,synthetic,setargs.K))
        data=[];
        
        for m in range(0,len(setargs.Docs)):
            
            data.append([m]);
        output = process_pool.starmap(parallel_assign, data)
        
        
        
        for m in range(0,len(setargs.Docs)):                    
                setargs.n_m_z[m]=output[m][0];
                setargs.n_m_t[m]=output[m][1]
                setargs.z_m_n.append(output[m][2])
                setargs.n_r[m]=output[m][3]
                setargs.doc_term[m]=output[m][4]
        process_pool.close()
        
        
        
        for i in range(setargs.V):
            setargs.n_z_t[:,i]=numpy.sum(numpy.sum(setargs.n_m_t[::,i::setargs.V,:],axis=1),axis=0)
        setargs.n_r=numpy.cumsum(setargs.n_r,axis=1);
        
    
def parallel_assign(m):
    
    
    doc=docs[m];
    n_m_z=_n_m_z[m].copy()
    n_m_t=_n_m_t[m].copy()
    n_r=_n_r[m]
    doc_term=_doc_term[m].copy();
    variants_dict=_variants_dict;
    z_n=[];
    t1=0
   
   
    
    #for na,t in enumerate(doc):   
    list_weight=[];
    ref_list=set();
    for r in doc:
                #if seed words in dictionary
                if(r in variants_dict): 
                   
                                    
                    
                    c1=[];
                    c1=variants_dict.get(r)
                    
                    
                    #print(c1,len(c1));
                    list1=[];
                    weightslist=[];
                    
                    p=100/(len(c1)+(K-known));
                    p1=(100-p)/len(c1);
                    weightslist=[p1] * len(c1);
                    list1.extend(c1);
                    if(K!=known):
                        weightslist.extend([p/(K-known)]*(K-known))
                        list1.extend(np.arange(known,K))
                    #z=np.random.choice(list1, size=1)[0]
                    #print(len(list1),len(weightslist))
                    z=random.choices(list1, weights=weightslist, k=1)[0]
                    
                   
                else:
                    z = numpy.random.randint(0, K)
                #z=numpy.random.multinomial(1, doc_topic[m]/ doc_topic[m].sum()).argmax()
                #print(z)

                doc_term[r]+=1;
                #print(setargs.doc_term[m,r])
                z_n.append(z)
                n_m_z[z] += 1
                n_m_t[r,z]+= 1
                #setargs.n_z_t[z, r] += 1
                #setargs.n_z[z] += 1
                t1+=1;
                rr=int(r/4)
                n_r[rr+1]+=1
                #setargs.n_r_max[m,rr+1]=max(n_r_max[m,rr],na)
                
                
    #z_m_n.append(numpy.array(z_n))

    return n_m_z,n_m_t,z_n,n_r,doc_term;   

    
def NMF_parallel_assign(m):
   
    doc=docs[m]
    n_m_z=_n_m_z[m]
    n_m_t=_n_m_t[m]
    n_r=_n_r[m]
    doc_term=_doc_term[m]
    theta_est=_theta_est[m]
    
    
    z_n=[];
    t1=0
    temp_array=np.zeros((V,K))
    for j in range(0,V):
            temp_array[j]=np.multiply(theta_est,word_topic[:,j])
    
    #for na,t in enumerate(doc):
       
    for r in doc:
                #if seed words in dictionary
            
                prob=temp_array[r];
                prob=numpy.asarray(prob)/numpy.sum(prob)
                prob=numpy.nan_to_num(prob, nan=0)
                z1 = numpy.nonzero(numpy.random.multinomial(1,prob,size = 1)[0])
                z=z1[0][0]  
                #z=numpy.random.multinomial(1, doc_topic[m]/ doc_topic[m].sum()).argmax()
                #print(z)

                doc_term[r]+=1;
                #print(setargs.doc_term[m,r])
                z_n.append(z)
                n_m_z[z] += 1
                n_m_t[r,z]+= 1
                #setargs.n_z_t[z, r] += 1
                #setargs.n_z[z] += 1
                t1+=1;
                rr=int(r/4)
   
  
                n_r[rr+1]+=1
                #setargs.n_r_max[m,rr+1]=max(n_r_max[m,rr],na)

    
    return n_m_z,n_m_t,z_n,n_r,doc_term;  


def NMF_parallel_initialize(variants_dict,H,W):
    
    #Initialization using NMF
    word_topic=numpy.zeros((setargs.K, setargs.V))
    theta_est=W/W.sum(axis=1)[:,None]

    #print(theta_est)

    for i in range(int(setargs.V/4)):

            n_z_t_sub=H[:,i*4:(i+1)*4]
            #print(n_z_t_sub)
            n_z_sub=numpy.transpose(numpy.sum(n_z_t_sub,axis=1))
            word_topic[:,i*4:(i+1)*4]=n_z_t_sub / n_z_sub[:, numpy.newaxis]
            word_topic[:,i*4:(i+1)*4]=numpy.nan_to_num(word_topic[:,i*4:(i+1)*4]) 
            #print(word_topic[:,i*4:(i+1)*4].shape)

    for i in range(setargs.V):
            word_topic[:,i]=word_topic[:,i]/word_topic[:,i].sum()

    word_topic=numpy.nan_to_num(word_topic)
    #print(word_topic)


    #pd.DataFrame(word_topic).to_csv(r'H.csv',index=False)
    #pd.DataFrame(theta_est).to_csv(r'W.csv',index=True)
    
    process_pool = multiprocessing.Pool(50,initializer = NMF_child_initialize, initargs = (setargs.Docs,variants_dict,setargs.n_m_z,setargs.n_m_t,setargs.n_r,setargs.doc_term,theta_est,word_topic,setargs.K,setargs.V))
    data=[];
   
    for m, doc in enumerate(setargs.Docs):
        data.append([m])
        
    print('start parallelization')
                            
    output = process_pool.starmap(NMF_parallel_assign, data)
    
    print('process output')
    for m in range(len(setargs.Docs)):
                setargs.n_m_z[m]=output[m][0];
                setargs.n_m_t[m]=output[m][1]
                setargs.z_m_n.append(output[m][2])
                setargs.n_r[m]=output[m][3]
                setargs.doc_term[m]=output[m][4]
        
    process_pool.close();
    
    for i in range(setargs.V):
            setargs.n_z_t[:,i]=numpy.sum(numpy.sum(setargs.n_m_t[::,i::setargs.V,:],axis=1),axis=0)
    setargs.n_r=numpy.cumsum(setargs.n_r,axis=1);
    
  
    filename = "Demixer_cpp/input/n_m_z.dat"
    fileobj = open(filename, mode='wb')
    setargs.n_m_z.tofile(fileobj)
    fileobj.close
    
    filename = "Demixer_cpp/input/n_z_t.dat"
    fileobj = open(filename, mode='wb')
    setargs.n_z_t.tofile(fileobj)
    fileobj.close
    
    filename = "Demixer_cpp/input/n_m_t.dat"
    fileobj = open(filename, mode='wb')
    setargs.n_m_t.tofile(fileobj)
    fileobj.close
    
