import numpy as np
import os
import pickle  
import seaborn as sns;
import pandas as pd;
import matplotlib.pyplot as plt
from plotnine import *
import plotnine as p9
import math
import pdb

def rela_error(x,y):
            x=x;
            y=y;
            z=abs(x-y)/(x+0.01);
            return np.mean(z);

rel_error_invitro=np.zeros((6,10))
ground_truth=pd.read_csv('../Figures/invitro/invitro_groundtruth.csv')
list3=np.array(ground_truth['major'].values)
list4=np.array(ground_truth['minor'].values)

path=["invitro_ref_Kset_nowt","invitro_ref_Kset_wtset","invitro_ref_highK_wtset","invitro_ref_highK_nowt","invitro_noref_highK","invitro_noref_Kset"]

for i in range(0,6):
    #os.system("python main_preprocess.py ~/LVM_Multistrain/Preprocessing/WGS_Indian_data/Malwai_samples/invitro/invitro_new_0.vcf.gz "+str(path[i])+" tbprof 2;")
    for j in range(0,10):
        
        os.system("./Demixer finaloutput/"+str(path[i])+"/")
        os.system("python postprocessing.py finaloutput/"+str(path[i])+"/")
        os.system("Rscript ../Figures/clustering.R ~/LVM_Multistrain/Preprocessing/Pycode/Demixer-main/Demixer/scripts/finaloutput/"+str(path[i])+"/")
        new_theta=pd.read_csv("finaloutput/"+str(path[i])+"/final_proportion.csv")
        
        final_max1=np.zeros(len(new_theta))
        final_max2=np.zeros(len(new_theta))

        for k in range(0,len(new_theta)):
            final_max1[k]=max(new_theta.iloc[k,:])
            final_max2[k]=min(new_theta.iloc[k,:])

        rel_error_invitro[i,j]=rela_error(list3,np.array(final_max1))+rela_error(list4,np.array(final_max2))/2
            
print(rel_error_invitro)

with open('ablation/invitro_ablation.pkl', 'wb') as f:
    pickle.dump(rel_error_invitro, f)
f.close();
