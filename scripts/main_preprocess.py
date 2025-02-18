#!/usr/bin/python
import re,sys
import gc
import os
import math
import time
import numpy
from pandas import ExcelWriter
from pandas import ExcelFile
import scipy.sparse as sps
import pickle
import pdb;
import random
from init import*
from vcf_preprocess import*
from scipy.stats import entropy
import re,sys
import pandas as pd;
import os
import numpy as np;
import random
from setargs import*
import psutil
import pickle


#Function to reset the variables
def reset():
    setargs.z_m_n=[];
    
    

    setargs.n_z_t = numpy.zeros((setargs.K, setargs.V),dtype='int32') # word count of each topic and vocabulary
    setargs.n_z = numpy.zeros(setargs.K)        # word count of each topic
    
    shared_n_m_t=multiprocessing.Array('i',len(setargs.Docs_1)*setargs.V*setargs.K, lock=False)
    shared_array = np.ctypeslib.as_array(shared_n_m_t)
    shared_array = shared_array.reshape(len(setargs.Docs_1),setargs.V,setargs.K)
    
    setargs.n_m_t=shared_array
    
    shared_n_m_z=multiprocessing.Array('i',len(setargs.Docs_1)*setargs.K, lock=False)
    shared_array1 = np.ctypeslib.as_array(shared_n_m_z)
    shared_array1 = shared_array1.reshape(len(setargs.Docs_1),setargs.K)
    
    setargs.n_m_z=shared_array1
    #setargs.n_m_z=0
    #setargs.n_m_t=0
    setargs.n_r=0
    setargs.doc_term=0
    


#Function to load the trained variables    
def load():
    with open('finaloutput/'+setargs.strain+'/n_m_z.pkl', 'rb') as f:
        setargs.n_m_z=pickle.load(f)
    f.close();
        
    with open('finaloutput/'+setargs.strain+'/n_z_t.pkl', 'rb') as f:
        setargs.n_z_t=pickle.load(f)
    f.close()
    
    with open('finaloutput/'+setargs.strain+'/z_m_n.pkl', 'rb') as f:
        setargs.z_m_n=pickle.load(f)
    f.close();
        
    with open('finaloutput/'+setargs.strain+'/n_m_t.pkl', 'rb') as f:
        setargs.n_m_t=pickle.load( f)
    f.close()
        
    with open('finaloutput/'+setargs.strain+'/n_r.pkl', 'rb') as f:
        setargs.n_r=pickle.load(f)

        
def construct_dict(doc,synthetic=False):
    variants_1=pd.read_csv('finaloutput/'+setargs.strain+'/subset.csv')
    col_name=setargs.Docs_col_name
    if synthetic==False:

            
            flat_list = [int(c) for c,xs in enumerate(doc) if xs!=0]
            
            v_row,v_col=variants_1.shape    
            variants_dict={}
            
            for i in range(v_col-1,-1,-1):
                l1=variants_1.iloc[:,i]  
                common_set=list(set(l1)&set(flat_list))
                
                if(len(common_set)>0):
                    for j in range(0,len(common_set)):
                                c1=common_set[j]
                                if c1 in variants_dict.keys():
                                    if(i not in variants_dict.get(int(c1))):
                                        variants_dict[int(c1)].append(i);
                                    #print(variants_dict.get(c1),count);
                                else:
                                    variants_dict.setdefault(c1, [])
                                    variants_dict[int(c1)].append(i);
                    variant_name=variants_1.columns[i];
                    for j in range(1,len(variant_name),2):
                        if(variant_name[0:j] in variants_1.columns):
                            l1=variants_1[variant_name[0:j]]
                            for k in range(0,len(l1)):
                                    c1=l1[k]
                                    if(not math.isnan(l1[k])):
                                        if c1 in variants_dict.keys():
                                            if(i not in variants_dict.get(int(c1))):
                                                variants_dict[int(c1)].append(i);
                                            for subv in range(j+2,len(variant_name),2):
                                                try:
                                                    if(list(variants_1.columns.values).index(variant_name[0:subv]) not in variants_dict.get(int(c1))):
                                                        variants_dict[int(c1)].append(list(variants_1.columns.values).index(variant_name[0:subv]))
                                                except ValueError :
                                                    u=0;


                                        else:
                                            variants_dict.setdefault(c1, [])
                                            variants_dict[int(c1)].append(i);   
                                            for subv in range(j+2,len(variant_name),2):
                                                try:
                                                    if(list(variants_1.columns.values).index(variant_name[0:subv]) not in variants_dict.get(int(c1))):
                                                        variants_dict[int(c1)].append(list(variants_1.columns.values).index(variant_name[0:subv]))
                                                except ValueError :
                                                    u=0;
    
    
    return variants_dict

def write_datfile(filename,fileobject):
    filepointer = open(filename, mode='wb')
    fileobject.tofile(filepointer)
    filepointer.close  

if __name__ == '__main__':
    
    #..............Preprocessing of vcf file........................................................

    st = time.time() 
    filename=sys.argv[1]
    print('reading start');
    if(int(sys.argv[4])==2):
        setargs.callset = allel.read_vcf(filename,fields=["variants/POS","calldata/AD","samples"],buffer_size=1048576,chunk_length=65536)
    elif(int(sys.argv[4])==3):
        setargs.callset = allel.read_vcf(filename,fields=["variants/POS","calldata/COV","samples"],buffer_size=1048576,chunk_length=65536)

    print('reading done');
    
    if(sys.argv[3]=='covid'):
        genotype=allel.vcf_to_dataframe(filename,fields=['POS','REF','ALT','is_snp']);
    else:
        genotype=cal_genotype(filename);
    
    
    
    genotype=genotype[genotype['is_snp']==True]
    genotype=genotype.reset_index();
    row,col=genotype.shape 
    print(row,col);
    
    read_countdata=pd.DataFrame()
    
    sample_name=list()
    length=len(setargs.callset['samples'])
    doc_n=length
    setargs.idf=[]
    for j in range(0,length):     
            sample_name.append(setargs.callset['samples'][j]);
  
    read_countdata['Sample_id']=sample_name;
    
    print("Number of samples: "+ str(length))
    setargs.Docs=[[] for x in range(length)]
    setargs.Docs_1=[[] for x in range(length)]
    setargs.Docs_2=[[] for x in range(length)]
    setargs.Docs_3=[[] for x in range(length)]
    setargs.Docs_col_name1=[];
    setargs.Docs_col_name2=[];
    setargs.Docs_col_name3=[];
    setargs.Docs_col_name=[]
    setargs.idf1=[]
    setargs.idf2=[]
    setargs.idf3=[]
    id_num=0
    weight=1
    
    
    #int(sys.argv[4])
    
    for variant in range(0,row): 
        read_data(variant,genotype,length,read_countdata,weight) 
    
    setargs.strain=sys.argv[2];    #output_filename
    current_directory = os.getcwd()
    out_dir = os.path.join(current_directory, r'finaloutput/'+setargs.strain)
    try:
        os.mkdir(out_dir)
    except OSError as error:
        print(error)
    
    pd.DataFrame(read_countdata).to_csv(r'finaloutput/'+setargs.strain+'/Sample_id.csv',index=False)
    setargs.Docs_col_name=setargs.Docs_col_name1+setargs.Docs_col_name2+setargs.Docs_col_name3
    setargs.idf=setargs.idf1+setargs.idf2+setargs.idf3
    pickle.dump(setargs.Docs_col_name, open('finaloutput/'+setargs.strain+'/Docs_col_name.pkl', 'wb'), protocol=4)
    
    
    setargs.V=len(setargs.Docs_col_name);   #vocabulary
    setargs.snptype=np.zeros(setargs.V,dtype='uint16')
    setargs.snppos=np.zeros(setargs.V,dtype='uint16')
    setargs.Docs=np.array(setargs.Docs_1)
    setargs.Docs_2=np.array(setargs.Docs_2)
    setargs.Docs_3=np.array(setargs.Docs_3)
    #pdb.set_trace();
    x1=setargs.Docs.shape[1]
    x2=setargs.Docs_2.shape[1]+setargs.Docs.shape[1]
    
    for j in range(0,setargs.V):
        if(j<setargs.Docs.shape[1]):
            setargs.snptype[j]=2;
            setargs.snppos[j]=0;
        elif(j<setargs.Docs_2.shape[1]+setargs.Docs.shape[1]):
            setargs.snptype[j]=3;
            setargs.snppos[j]=setargs.Docs.shape[1];
        else:
            setargs.snptype[j]=4;
            setargs.snppos[j]=setargs.Docs_2.shape[1];
    #pdb.set_trace();
    setargs.Docs_1=np.concatenate((setargs.Docs, setargs.Docs_2), axis=1)
    setargs.Docs_1=np.concatenate((setargs.Docs_1, setargs.Docs_3), axis=1)
    print(setargs.Docs.shape[1],setargs.Docs_2.shape[1]+setargs.Docs.shape[1],setargs.Docs_1.shape[1]); 	   
    #pdb.set_trace();

    setargs.Docs_1=setargs.Docs_1.astype("int16")
    write_datfile("finaloutput/"+setargs.strain+"/Docs_int_onlySnpgap.dat",setargs.Docs_1)
    setargs.V=len(setargs.Docs_col_name);   #vocabulary
    setargs.cols_name=dict(zip(setargs.Docs_col_name,range(len(setargs.Docs_col_name))))
    

    variants_dict={}
     
    variants=pd.read_csv('db/'+sys.argv[3]+'_snpset.csv')        #for invitro
    
    
    #variants=pd.DataFrame()   #included for ablation study - to be removed
    
    (v_row, v_col)=variants.shape
    
    for grid_no in range(1,2):

        mut_l=[];
        
        #Addig seed words on to dictionary: Only those strains whose seed words are present in the dataset 
        #would be considered for further analysis
        sum=0;
        print('variants dictionary')

        variants_1=pd.read_csv('db/'+sys.argv[3]+'_unique.csv')  
      
        
        #variants_1=pd.DataFrame()   #included for ablation study - to be removed
        
        variants_subset=pd.DataFrame();
        variant_ids_subset=pd.DataFrame();
        variants_subset1=pd.DataFrame();
        variant_ids_subset1=pd.DataFrame();
        count=0;
        col_names=[]
        for i in range(0,v_col):
            if(1):
            
            #if(variants_1.columns[i]!='bov'):
                l1=variants_1.iloc[:,i]
                sum=sum+len(list(set(l1)&set(setargs.cols_name)))
                sub_list=[]
                sub_list1=[];
                sub_ids_list=[]
                sub_ids_list1=[];
                l2=variants.iloc[:,i]   #to be removed
               
                if(len(list(set(l1)&set(setargs.cols_name)))>0):
                    #setargs.alpha.append(Possible_alpha[i])
                    print(variants_1.columns[i])
                    for j in range(0,len(l1)):
                        if(l1[j] in setargs.cols_name):
                            sub_list.append(setargs.cols_name.get(l1[j]))
                            sub_ids_list.append(l1[j])

                    for j in range(0,len(l2)):
                        if(l2[j] in setargs.cols_name):
                            sub_list1.append(setargs.cols_name.get(l2[j]))
                            sub_ids_list1.append(l2[j])

                    if((len(l1)>=5 and len(sub_list)>=5)):      

                        variants_subset = pd.concat([variants_subset,pd.Series(sub_list)], ignore_index=True, axis=1)
                        variants_subset1 = pd.concat([variants_subset1,pd.Series(sub_list1)], ignore_index=True, axis=1)
                        variant_ids_subset = pd.concat([variant_ids_subset,pd.Series(sub_ids_list)], ignore_index=True, axis=1)
                        variant_ids_subset1 = pd.concat([variant_ids_subset1,pd.Series(sub_ids_list1)], ignore_index=True, axis=1)
                        
                        col_names.append(variants_1.columns[i])
                        l1=variants.iloc[:,i]
                        for j in range(0,len(l1)):
                                
                                if(l1[j] in setargs.cols_name):
                                    c1=setargs.cols_name.get(l1[j])
                                    if c1 in variants_dict.keys():
                                        variants_dict[c1].append(count);
                                        #print(variants_dict.get(c1),count);
                                    else:
                                        variants_dict.setdefault(c1, [])
                                        variants_dict[c1].append(count);

                        count=count+1;
                        
        #pdb.set_trace()
        variants_subset.columns= col_names  
        variants_subset1.columns= col_names  
        variant_ids_subset.columns= col_names  
        variant_ids_subset1.columns= col_names  
        #print(variants_dict)            
        print(sum,count);
        
        pd.DataFrame(variants_subset).to_csv(r'finaloutput/'+setargs.strain+'/subset.csv',index=False) 
        
        pd.DataFrame(variants_subset1).to_csv(r'finaloutput/'+setargs.strain+'/subset1.csv',index=False) 
        pd.DataFrame(variant_ids_subset).to_csv(r'finaloutput/'+setargs.strain+'/subset_ids.csv',index=False) 
        
        pd.DataFrame(variant_ids_subset1).to_csv(r'finaloutput/'+setargs.strain+'/subset1_ids.csv',index=False) 

        f=open('finaloutput/'+setargs.strain+'/temp1.txt', "w")
        for j, key in enumerate(variants_dict):
            l=variants_dict.get(key)
            s1=""
            if(len(l)>0):               #if(len(l)>0):
                s1 = " ".join(str(e) for e in l)
            f.write(str(int(key))+ ' '+str(len(l))+' ' +s1+'\n')
        f.close()
        

        print(setargs.V,'max',np.max(setargs.Docs_1))
       
        #############ablation study need to be uncommented and commented###########
       
        setargs.K=count+2;  
        #setargs.K=30       #invitro K - higher value
        #setargs.K=19       #invitro K - set
        ##########################################################################

        with open('finaloutput/'+setargs.strain+'/'+'config.txt', 'w') as file:
            file.write(str(len(setargs.Docs_1)) + ' '+str(setargs.K)+' '+str(setargs.V)+' '+str(x1)+' '+str(x2));
        
        
        print('Initialization started')
        #pdb.set_trace()
        print(len(sys.argv))
        if(len(sys.argv)==6):
            weight=parallel_initialize(variants_dict,count,False,int(sys.argv[5]));     
        else:
            if(sys.argv[3]=='quanttb'):
                weight=parallel_initialize(variants_dict,count,True);               #quanttb
            else:
                weight=parallel_initialize(variants_dict,count,False);              #tbprof
            
        
        print(weight) 

    
        setargs.idf=np.array(setargs.idf)
        setargs.idf=setargs.idf.astype("int16")
        
       
        
        ############ manuscript revision ablation study##############
        setargs.idf[setargs.idf>0]=weight           #uncomment otherwise
        #setargs.idf[setargs.idf>0]=0
        ###############################################
        
        
        setargs.alpha=np.full(setargs.K,1,dtype='<f4')
        setargs.beta=np.full(setargs.K,0.01,dtype='<f4')
        write_datfile("finaloutput/"+setargs.strain+"/idf.dat",setargs.idf)
        write_datfile("finaloutput/"+setargs.strain+"/alpha.dat",setargs.alpha)
        write_datfile("finaloutput/"+setargs.strain+"/beta.dat",setargs.beta)
        print("time",time.time()-st)
        process = psutil.Process()
        memory_info = process.memory_info()  
        memory_used_gb = memory_info.rss / (1024 ** 3)
        print(f"Memory used: {memory_used_gb:.2f} GB")          
      
        '''
        #######parallelcheck
        
        cores=[1,2,4,8,16,32,64]
        timelist=[]
        for i in cores:
            st=time.time()
            os.system("./ablation/Demixer_parallel1 finaloutput/parallel1 "+str(i))    
            timelist.append(time.time()-st)
        with open('invitro_cores_list1.pkl', 'wb') as file:
            pickle.dump(timelist, file)    
        print(timelist)
        

        
        timelist=[]
        for i in cores:
            st=time.time()
            os.system("./ablation/Demixer_parallel4 finaloutput/parallel4 "+str(i))    
            timelist.append(time.time()-st)
        with open('invitro_cores_list4.pkl', 'wb') as file:
            pickle.dump(timelist, file)    
        print(timelist)
            
        timelist=[]
        for i in cores:
            st=time.time()
            os.system("./ablation/Demixer_parallel5 finaloutput/parallel5 "+str(i))    
            timelist.append(time.time()-st)
        with open('invitro_cores_list5.pkl', 'wb') as file:
            pickle.dump(timelist, file)    
        print(timelist)
       
        '''
        