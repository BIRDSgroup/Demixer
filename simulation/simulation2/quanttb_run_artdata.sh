
#!/bin/bash

#${1} - type {tbprof,quanttb}
#{2} - coverage
if [ "${1}" == "tbprof" ]; then
    ref_name="tuberculosis.fna"
else
    ref_name="CP003248.2.fasta"
fi

echo ${1} ${2}
do_one() { pilon --genome Pycode/Art_simulator/ref/${3} --bam Pycode/Art_simulator/data_${1}_coverage_${2}/${4}.sorted.bam --output art_${4} --outdir Pycode/Art_simulator/data_${1}_coverage_${2}/temp --vcf  --threads 40 --fix none; python quanttb_copy/quanttb/scripts/quanttb.py quant -v Pycode/Art_simulator/data_${1}_coverage_${2}/temp/art_${4}.vcf -o Pycode/Art_simulator/data_${1}_coverage_${2}/quanttb_output/${4}.txt; }

export -f do_one;

do_tbprof() { pilon --genome Pycode/Art_simulator/ref/${3} --bam Pycode/Art_simulator/data_${1}_coverage_${2}/${4}.sorted.bam --output art_${4} --outdir Pycode/Art_simulator/data_${1}_coverage_${2}/temp --vcf  --threads 40 --fix none; python quanttb_copy/quanttb/scripts/quanttb.py quant -v Pycode/Art_simulator/data_${1}_coverage_${2}/temp/art_${4}.vcf -db Pycode/Art_simulator/ref_tbprof/qtb.txt.db -o Pycode/Art_simulator/data_${1}_coverage_${2}/quanttb_output/${4}.txt; }

export -f do_tbprof

mkdir Pycode/Art_simulator/data_${1}_coverage_${2}/temp
mkdir Pycode/Art_simulator/data_${1}_coverage_${2}/quanttb_output

if [ "${1}" == "tbprof" ]; then
    seq 1 200 | parallel -P 5 do_tbprof ${1} ${2} ${ref_name} 
else
    seq 1 200 | parallel -P 5 do_one ${1} ${2} ${ref_name} 
fi

