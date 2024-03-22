import numpy as np
import pandas as pd
import csv
import numpy as np
import math
from scipy import stats
global K
global V
global Docs
global Docs_col_name
global cols_name
global q;
global strain
global z_m_n
global variants_dict
global n_m_z
global n_z_t
global n_z# word count of each topic
global n_m_t
global n_r
global n_r_max
global doc_term
global alpha
global beta
global z_n
global word;
global sampledata;
global variant_dict;
global callset;
global sampledata;
global idf;
global Docs_1;
global snptype
global snppos
global seeds
global Docs_test

def write_datfile(filename,fileobject):
    filepointer = open(filename, mode='wb')
    fileobject.tofile(filepointer)
    filepointer.close  

def rela_error(x,y):
            x=x+1;
            y=y+1;
            
            z=abs(((x-y))/(x));
            #print(z)
            #print(np.mean(z))
            return np.mean(z);
        
def rela_error(x,y): 
            z=abs(((x-y))/(x+0.01));
            return np.mean(z);
        
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def plot_prop_major(data,snp1,snp2,l):
    readcount=data
    total=0
    step=2;
    prop1=[]
    prop2=[]
    snp_num1=[]
    snp_num2=[]
    count=1;
    n_count=0;
    for i in range(0,snp1,step):
        snp_num1.append(count)
        snp_num2.append(count)
        count+=1;
        total=total+np.sum(readcount[i:i+step])
        n_count+=1;
        if(1):
            a=readcount[i]/total
            b=np.sum(readcount[i+1:i+step])/total
       
            if(a<b):
                if(a!=0 and (not math.isnan(a))):
                    #print(Docs_col_name[i])
                    prop1.append(a)
            else:
                if(b!=0 and (not math.isnan(b))):
                    #print(Docs_col_name[i])
                    prop1.append(b)
        total=0
    step=3;
    for i in range(snp1,snp2,step):
        snp_num1.append(count)
        snp_num2.append(count)
        count+=1;
        total=total+np.sum(readcount[i:i+step])
        #if(readcount[i]/total>0.05 and readcount[i]/total<0.95):
        if(1):
            a=readcount[i]/total
            b=np.sum(readcount[i+1:i+step])/total
            if(a<b):
                if(a!=0 and (not math.isnan(a))):
                    prop1.append(a)
            else:
                if(b!=0 and (not math.isnan(b))):
                    prop1.append(b)
        total=0
        
    step=4;
    for i in range(snp2,l,step):
        snp_num1.append(count)
        snp_num2.append(count)
        count+=1;
        total=total+np.sum(readcount[i:i+step])
        #if(readcount[i]/total>0.05 and readcount[i]/total<0.95):
        if(1):
            a=readcount[i]/total
            b=np.sum(readcount[i+1:i+step])/total
            if(a<b):
                if(a!=0 and (not math.isnan(a))):
                    prop1.append(a)
            else:
                if(b!=0 and (not math.isnan(b))):
                    prop1.append(b)
        total=0
        
    #print(prop1[prop1>0.03 and prop1<0.92])
    #plt.scatter(snp_num1+snp_num2, prop1+prop2, color='blue', label='prop1')
    #plt.hist(prop1,bins=10)
    
    data, bins = np.histogram(prop1, bins=10)
    max_count_bin = np.argmax(data)

    max_count = data[max_count_bin]
    
    #plt.figure(figsize=(16, 10))
    
    # Plot histogram for reference
    #histogram = sns.histplot(prop1, bins=10, )
    mode=stats.mode(np.around(prop1,2))
    #data, bins, _ = plt.hist(prop1, bins=10,)
    
    bin_index = np.argsort(data)[-1:][::-1]
    
    freq=[]
    if(data[bin_index]>0):
        if(bin_index>0):
            freq.append(bins[bin_index-1])

        if(bin_index<len(data)-1):
            freq.append(bins[bin_index+1])
        freq.append(bins[bin_index])
        
        m=np.mean(freq)
    else:
        m= 0
       
    if len(mode[0])>0:
        return mode[0][0],max_count,m
    else:
        return 0,0,0
    
def write_csvfile(data):
    max_length = max(len(row) for row in data)
    for row in data:
        while len(row) < max_length:
            row.append(None)
    return data

def thresh_conf(dir,test=False):
    Demixer_result=[]     
    
    if(test==True):
        diverg=list(pd.read_csv('finaloutput/Cryptic_rerun/diverg_id.csv')['0'].values)
        snp_conf=pd.read_csv('finaloutput/Cryptic_rerun/snp_plot_conf.csv')
    else:
        snp_conf=pd.read_csv(dir+'snp_plot_conf.csv')
        diverg=list(pd.read_csv(dir+'diverg_id.csv')['0'].values)
    data1=snp_conf.long_bar.values
    
    percentiles = [25, 75]
    col=['red','green']
    thresholds = np.percentile(data1, percentiles)
    prop_data=pd.read_csv(dir+'Prior.csv')

    with open(dir+'prop.csv', 'w', newline='') as file1, open(dir+'lineage.csv', 'w', newline='') as file2, open(dir+'Mixpred.csv', 'w', newline='') as file3:
        writer1 = csv.writer(file1)
        writer2 = csv.writer(file2)
        writer3 = csv.writer(file3)
        final_p=[]
        final_w=[]
        conf=[]
        for i in range(0,len(prop_data)):
            
            prop_list=prop_data.iloc[i,:][prop_data.iloc[i,:]>0].values
            id_list=np.array(diverg)[np.where(prop_data.iloc[i,:]>0)]
            words_list=[]
            for j in id_list:
                #words=j.split("/")
                words_list.append(j)

            if(len(words_list)>1):  

                            theta_prop=[0]*len(prop_list)

                            strain_ids=words_list.copy()
                            strain_ids.sort()

                            ran=0;
                            while(ran<len(strain_ids)):
                                if(words_list.count(strain_ids[ran])==1):
                                    ind=words_list.index(strain_ids[ran])
                                    theta_prop[ran]=prop_list[ind]
                                    ran+=1;
                                else: 
                                    ind=[index for index, value in enumerate(words_list) if value == strain_ids[ran]]
                                    for ran2 in ind:
                                        theta_prop[ran]=prop_list[ran2]
                                        ran+=1


                            final_prop=[]
                            final_strains=[]

                            for j in range(0,len(strain_ids)):
                                if(j!=len(strain_ids)-1):
                                    flag=0;

                                    if(strain_ids[j+1].startswith(strain_ids[j])):

                                        flag=1;
                                        if(theta_prop[j+1]>theta_prop[j]):
                                                strain_ids[j+1]=strain_ids[j+1]
                                        else:
                                                strain_ids[j+1]=strain_ids[j]
                                        theta_prop[j+1]=theta_prop[j]+theta_prop[j+1]

                                    else:


                                        final_strains.append(strain_ids[j])
                                        final_prop.append(theta_prop[j])    
                                        if(j+1==len(strain_ids)-1):
                                            final_strains.append(strain_ids[j+1])
                                            final_prop.append(theta_prop[j+1])  

                                    if(flag==1 and j+1==len(strain_ids)-1):
                                        final_strains.append(strain_ids[j+1])
                                        final_prop.append(theta_prop[j+1])

                                if(len(strain_ids)==1):
                                    final_strains.append(strain_ids[j])
                                    final_prop.append(theta_prop[j])

                            if(len(final_prop)==1):

                                Demixer_result.append(["Non-Mix"])

                            else:
                                if('_' in final_strains):

                                    if(snp_conf['long_bar'][i]>thresholds[1] ):          #90 for cryptic, 250 for Malawi

                                        Demixer_result.append(["Mix"])   
                                    else:
                                        final_strains_copy=[x for x in final_strains if x !='_']

                                        if(len(final_strains_copy)==1):
                                            Demixer_result.append(["Non-Mix"])
                                            final_prop=[1]
                                            final_strains=final_strains_copy
                                        else:
                                            Demixer_result.append(["Mix"])
                                          
                                else:
                                    Demixer_result.append(["Mix"]) 
                            final_w.append(final_strains)
                            final_p.append(final_prop)

            else:
                final_w.append(words_list)
                final_p.append(list(prop_list))

                Demixer_result.append(["Non-Mix"])

            if(Demixer_result[-1][0]=="Mix"):
                if(snp_conf['long_bar'][i]<thresholds[0]):
                    conf.append("low")
                elif(snp_conf['long_bar'][i]>=thresholds[0] and snp_conf['long_bar'][i]<=thresholds[1]):
                    conf.append("medium")
                elif(snp_conf['long_bar'][i]>thresholds[1]):
                    conf.append("high")
            #else:
            #    conf.append("")

            else:
                if(snp_conf['long_bar'][i]<thresholds[0]):
                    conf.append("high")
                elif(snp_conf['long_bar'][i]>=thresholds[0] and snp_conf['long_bar'][i]<=thresholds[1]):
                    conf.append("medium")
                elif(snp_conf['long_bar'][i]>thresholds[1]):
                    conf.append("low")
            #print(snp_conf['long_bar'][i],Demixer_result[-1],conf[-1])
            #pdb.set_trace()

            print(final_w[-1],final_p[-1],conf[-1])
        final_p=write_csvfile(final_p)
        final_w=write_csvfile(final_w)
       
        writer1.writerows(final_p)
        writer2.writerows(final_w)
        writer3.writerows(write_csvfile(Demixer_result))
        pd.DataFrame(conf).to_csv(dir+'confidence_new.csv',index='False')


