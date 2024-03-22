#!/bin/bash

#ref: Splitstrains code
set -e

#ref_name="tuberculosis.fna"

if [ "${3}" == "tbprof" ]; then
    ref_name="tuberculosis.fna"
else
    ref_name="CP003248.2.fasta"
fi
ref_folder="ref_"${3}
no_strain=$4
o_folder="data_"${3}"_coverage_"$2


#get the strain names

strain_names=$(head -n ${1} ${3}_${no_strain}_${2}_strain.txt | tail -n 1)

echo $strain_names

echo "generate reads"

touch $1_temp_1
touch $1_temp_2
var=${1}' '
if [ "${2}" == "70_30" ]; then
    arr=(70 30)
elif [ "${2}" == "90_10" ]; then
    arr=(90 10)
fi

echo ${arr[$((0))]}


for i in $(seq 1 $no_strain); do
    strain=$(echo $strain_names | cut -d ' ' -f $i);
    
    var=$var$(echo ${strain%.fasta*})' '
    
    
    #printf  ${var}"\t" >>strain_names.txt
    
    # -f $2   # -f ${arr[$((i-1))]}
    
    if [[ "$2" != *_* ]]; then
        
        art_illumina -ss HS25 -sam -i $ref_folder/$strain -l 101 -o $o_folder/strain_$1_${i}_R -qL 21 -qU 31 -qs -9 -qs2 -8 -na -f $2 -p -s 10 -m 300  
    else
        art_illumina -ss HS25 -sam -i $ref_folder/$strain -l 101 -o $o_folder/strain_$1_${i}_R -qL 21 -qU 31 -qs -9 -qs2 -8 -na -f ${arr[$((i-1))]} -p -s 10 -m 300  
    fi
   
    cat $1_temp_1  $o_folder/strain_$1_${i}_R1.fq >$1_temp
    mv $1_temp $1_temp_1
    cat $1_temp_2  $o_folder/strain_$1_${i}_R2.fq >$1_temp
    mv $1_temp $1_temp_2
    
    rm $o_folder/strain_$1_${i}_R1.fq 
    rm $o_folder/strain_$1_${i}_R2.fq
    rm $o_folder/strain_$1_${i}_R.sam
done
echo $var >>strain_names_${3}_coverage_$2.txt
mv $1_temp_1 $o_folder/$1_R1.fq
mv $1_temp_2 $o_folder/$1_R2.fq


echo "generate sam" $1

bwa mem -K 100000000 -Y -R  '@RG\tID:sample_'$1'\tLB:sample_'$1'\tPL:ILLUMINA\tPM:HISEQ\tSM:sample_'$1 ref/$ref_name $o_folder/$1_R1.fq $o_folder/$1_R2.fq > $o_folder/$1.aligned_reads.sam; 

rm $o_folder/$1_R1.fq

rm $o_folder/$1_R2.fq
# Sort the BAM file
samtools sort $o_folder/$1.aligned_reads.sam -o $o_folder/$1.sorted.bam

# Index the sorted BAM file
samtools index $o_folder/$1.sorted.bam

rm $o_folder/$1.aligned_reads.sam

freebayes -f ref/$ref_name  $o_folder/$1.sorted.bam > $o_folder/art_$2_$1.vcf
bgzip $o_folder/art_$2_$1.vcf
bcftools index $o_folder/art_$2_$1.vcf.gz

