import numpy as np
import os
import pickle  
import seaborn as sns;
import pandas as pd;
import matplotlib.pyplot as plt
from plotnine import *
import plotnine as p9
import math
from setargs import*
import pdb


rel_error_invitro=np.zeros((11,100))
ground_truth=pd.read_csv('../Figures/invitro/invitro_groundtruth.csv')
list3=np.array(ground_truth['major'].values)
list4=np.array(ground_truth['minor'].values)


for i in range(0,11):
    os.system("python main_preprocess.py ~/LVM_Multistrain/Preprocessing/WGS_Indian_data/Malwai_samples/invitro/invitro_new_0.vcf invitro_new_const"+str(i)+" tbprof 2 "+str(i))
    for j in range(0,100):
        os.system("./Demixer finaloutput/invitro_new_const"+str(i))
        os.system("python postprocessing.py finaloutput/invitro_new_const"+str(i)+"/")
        os.system("Rscript ../clustering.R ~/LVM_Multistrain/Preprocessing/Pycode/Demixer-main/Demixer-main/scripts/finaloutput/invitro_new_const"+str(i)+"/")
        new_theta=pd.read_csv("../scripts/finaloutput/invitro_new_const"+str(i)+"/final_proportion.csv")
        final_max1=np.zeros(len(new_theta))
        final_max2=np.zeros(len(new_theta))

        for k in range(0,len(new_theta)):
            final_max1[k]=max(new_theta.iloc[k,:])
            final_max2[k]=min(new_theta.iloc[k,:])

        rel_error_invitro[i,j]=rela_error(list3,np.array(final_max1))+rela_error(list4,np.array(final_max2))/2
pdb.set_trace()
with open('invitro.pkl', 'wb') as f:
    pickle.dump(rel_error_invitro, f)
f.close();