import numpy as np
import pandas as pd
import sys
import pdb
import pickle
from setargs import *
dir=sys.argv[1]
file_path = "finaloutput/Cryptic_rerun/config.txt"
filename = dir+"Docs_test.dat"
fileobj = open(filename, mode='r')
Docs_1 = np.fromfile(fileobj, dtype=np.int16)
fileobj.close
with open('finaloutput/Cryptic_rerun/Docs_col_name.pkl', 'rb') as f:
    Docs_col_name=pickle.load(f)
f.close();

with open(file_path, "r") as file:
    line = file.readline()
    integer_strings = line.split()
    integers = [int(i) for i in integer_strings]

K=integers[1]
V=integers[2]
snp1=integers[3]
snp2=integers[4]

with open(dir+"testconfig.txt", "r") as file:
    line = file.readline()
    integer_strings = line.split()
    integers = [int(i) for i in integer_strings]

m=integers[0]
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


pd.DataFrame(theta).to_csv(dir+'Prior.csv',index=False)


data=pd.DataFrame(Docs_1,columns=Docs_col_name)
cols_name=dict(zip(Docs_col_name,range(len(Docs_col_name)))) 

thresh_conf(dir, True)