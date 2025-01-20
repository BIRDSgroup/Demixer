import setargs;
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
import pdb
import pickle

#multiprocessing.set_start_method('spawn')
variant_dict_array=[]
r=0
#warnings.filterwarnings("ignore")
#os.environ['OPENBLAS_NUM_THREADS'] = '1'

weight=0;    
def output_to_file(dict_list,foldername,D):
    f = open(foldername+'temp2.txt', "w")
    #f = open("Snpgaponlysubvariant.txt", "w")
    total=0;
    count=0
    for i,item in enumerate(dict_list):

        if(item is not None):
            total=total+len(item)
            #print(len(item))
            count+=1
            f.write(str(i)+' '+str(len(item))+'\n');
            for j, key in enumerate(item):
                l=item.get(key)
                s1=""
                if(len(l)>0):
                    s1 = " ".join(str(e) for e in l)
                f.write(str(int(key))+ " "+str(len(l))+" "+s1+'\n')
    f.close()
    
    avg=float(total)/count 
    
    if(D!=0):
        #a=(setargs.Docs_1.shape[1]/4-avg)/setargs.K+0.0001 #for synthetic dataset
        a=(setargs.Docs_1.shape[1]/2-avg)/setargs.K+0.0001
        b=avg/D+0.0001
        weight=np.floor(a/b);
        print("weight",D, setargs.Docs_1.shape[1],setargs.Docs_1.shape[1]/4,avg,a,b,weight);
        if(a<b):
            weight=1
    else:
        weight=0

    weight=0
    
        
    print(weight)
    #pdb.set_trace()
    print("Average",total/count)
    
    return weight;
    

def parallel_construct(synthetic,foldername,fl,D,variants=pd.DataFrame()):

    process_pool = multiprocessing.Pool(100)
    data=[]
    
    dict_array=[]
    for m in range(0,len(setargs.Docs_1)):

            if(fl==1):
                data.append([m,synthetic,fl]);
            else:
                data.append([m,synthetic,fl,variants]);
            #dict_array.append(construct_dict(m,setargs.Docs_col_name,variants_1))
    out=process_pool.starmap(construct_dict, data)
    for m in range(len(setargs.Docs_1)):
        dict_array.append(out[m])

    weight=output_to_file(dict_array,foldername,D)

    process_pool.close()
    
    v=[dict_array,weight]
    return(v)
    #:process_pool.close()
    
def construct_dict(m,synthetic,fl,variants=pd.DataFrame()):
    
    doc=list(setargs.Docs_1[m])
    col_name=setargs.Docs_col_name
    if(fl==1):
        variants_1=pd.read_csv('finaloutput/'+setargs.strain+'/subset.csv')
        variants_2=pd.read_csv('finaloutput/'+setargs.strain+'/subset1.csv')
        s_name=dict(zip(variants_2.columns.values,range(len(variants_2.columns))))
    else:
        variants_1=variants
    flag=0;
    if synthetic==False:
          
            flat_list = [int(c) for c,xs in enumerate(doc) if xs!=0]
            #flat_list=doc
            
            v_row,v_col=variants_1.shape    
            variants_dict={}
            
            #variants_dict=Dict.empty(key_type=types.int32, value_type=nb.typed.List.empty_list(nb.types.int8),)
            
            for i in range(v_col-1,-1,-1):
                l1=variants_1.iloc[:,i]  
                common_set=list(set(l1)&set(flat_list))
                
                if(len(common_set)>0):
                    flag+=1;
                    for j in range(0,len(common_set)):
                                c1=common_set[j]
                                if c1 in variants_dict.keys():
                                    if(i not in variants_dict.get(int(c1))):
                                        variants_dict[int(c1)].append(i);
                                    #print(variants_dict.get(c1),count);
                                else:
                                    variants_dict.setdefault(c1, [])
                                    variants_dict[int(c1)].append(i);

                    variant_name=str(variants_1.columns[i]);
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
                                            #if(i not in variants_dict.get(int(c1))):
                                            variants_dict[int(c1)].append(i);   
                                            for subv in range(j+2,len(variant_name),2):
                                                try:
                                                    if(list(variants_1.columns.values).index(variant_name[0:subv]) not in variants_dict.get(int(c1))):
                                                        variants_dict[int(c1)].append(list(variants_1.columns.values).index(variant_name[0:subv]))
                                                except ValueError :
                                                    u=0;
    else:
            flat_list = [int(c) for c,xs in enumerate(doc) if xs!=0]
            v_row,v_col=variants_1.shape    
            variants_dict={}
            
            #variants_dict=Dict.empty(key_type=types.int32, value_type=nb.typed.List.empty_list(nb.types.int8),)
            
            for i in range(v_col-1,-1,-1):
                
                l1=variants_1.iloc[:,i]  
                common_set=list(set(l1)&set(flat_list))
                
                if(len(common_set)>0):
                    flag+=1;
                    for j in range(0,len(common_set)):
                                c1=common_set[j]
                                if c1 in variants_dict.keys():
                                    if(i not in variants_dict.get(int(c1))):
                                        variants_dict[int(c1)].append(i);
                                    #print(variants_dict.get(c1),count);
                                else:
                                    variants_dict.setdefault(c1, [])
                                    variants_dict[int(c1)].append(i);
                    
                    l2=variants_2.iloc[:,i]  
                    common_set=list(set(l2)&set(flat_list))

                    for j in range(0,len(common_set)):
                        c1=common_set[j]
                        if c1 in variants_dict.keys():
                            if(i not in variants_dict.get(int(c1))):
                                variants_dict[int(c1)].append(i);
                            #print(variants_dict.get(c1),count);
                        else:
                            variants_dict.setdefault(c1, [])
                            variants_dict[int(c1)].append(i);
                    
    #if(flag==0):
    #    print('no_common ',m)
        
    return variants_dict;

def parallel_initialize(variants_dict,known,synthetic=False,D=3):
        
        global variant_dict_array  
        v=parallel_construct(synthetic,'finaloutput/'+setargs.strain+'/',1,D)
        variant_dict_array=v[0]
        weight=v[1]
        return weight;
        
