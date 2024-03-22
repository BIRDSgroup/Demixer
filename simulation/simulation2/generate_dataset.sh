#!/bin/bash

#${1} - type {tbprof,quanttb}
#${2} - number pf strains
#${3} - coverage

type=${1}
coverage=${3}
if [ "${type}" == "tbprof" ]; then
    ref_name="tuberculosis.fna"
else
    ref_name="CP003248.2.fasta"
fi

do_one() { 
echo ${1} ${2} ${3} ${4} 
freebayes -f ref/${2} -@ art_${3}_quanttb.vcf.gz ${4}.sorted.bam > art_${3}_${4}.vcf;  
rm art_${3}_${4}.vcf.gz;
rm art_${3}_${4}.vcf.gz.csi;
bgzip art_${3}_${4}.vcf; 
bcftools index art_${3}_${4}.vcf.gz;
}

rm ${1}_${2}_${3}_strain.txt
python find_strain_names.py ${2} ${1}>> ${1}_${2}_${3}_strain.txt

mkdir data_${1}_coverage_${3}

echo ${3}"_coverage"

seq 1 200 |parallel -P 10 ./gen_vcf_nonparallel.sh {} ${3} ${1} ${2}

cd data_${1}_coverage_${3}


bcftools merge *.vcf.gz -Oz -o art_${3}_quanttb.vcf
bgzip art_${3}_quanttb.vcf
bcftools index art_${3}_quanttb.vcf.gz

export -f do_one
seq 1 200 | parallel -P 10 do_one ${type} ${ref_name} ${coverage} 

bcftools merge art_${3}_[0-9]*.vcf.gz -Oz -o art_${3}_quanttb.vcf
