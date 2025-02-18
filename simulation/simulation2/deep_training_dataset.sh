#generate art simulator generated reads for traing DL methods

#!/bin/bash

#ref: Splitstrains code
set -e

ref_name="tuberculosis.fna"

ref_folder="ref_tbprof"
o_folder="deep_train_sub_new"

echo "generate reads"

my_array=("1" "2" "2.1" "2.2.2" "3" "3.1" "3.1.2" "3.1.1" "3.1.2.1" "3.1.2.2")
for s in ${my_array[@]}; do
        strain="ref_tbprof/lineage$s.fasta"; 
        echo $strain
        id=$(echo $strain | sed 's/\.fasta$//'| cut -d "/" -f 2)
        echo $id
        art_illumina -ss HS25 -sam -i $strain -l 150 -p -o $o_folder/${id}_R -na -f 50 -m 200 -s 10 -s 0.5 -ef;
done





