import sys
from vcf_preprocess import*
import pickle
import pdb
import os

filename=sys.argv[1]
setargs.callset = allel.read_vcf(filename,fields=["variants/POS","calldata/AD","samples"],buffer_size=1048576,chunk_length=65536)


with open('finaloutput/Cryptic_rerun/Docs_col_name.pkl', 'rb') as f:
    setargs.Docs_col_name=pickle.load(f)
f.close();

setargs.cols_name=dict(zip(setargs.Docs_col_name,range(len(setargs.Docs_col_name))))

genotype=cal_genotype(filename);
genotype=genotype[genotype['is_snp']==True]
genotype=genotype.reset_index();
row,col=genotype.shape 
length=len(setargs.callset['samples'])
setargs.Docs_test=np.zeros((length,148990),dtype="int16")
for variant in range(0,row): 
        read_data_test(variant,genotype,length,[],0) 
        
out_dir = 'finaloutput/test'
try:
    os.mkdir(out_dir)
except OSError as error:
    print(error)      

filename = "finaloutput/test/Docs_test.dat"
fileobj = open(filename, mode='wb')
setargs.Docs_test.tofile(fileobj)
fileobj.close
        
with open('finaloutput/test/'+'testconfig.txt', 'w') as file:
    file.write(str(len(setargs.Docs_test)));
#pdb.set_trace()
print(setargs.Docs_test)
