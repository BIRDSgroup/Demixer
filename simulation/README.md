## Simulation

This folder has the scripts used for generating the different synthetic data and running Demixer on them.

### Simulation 1
Run ```python simulate_dataset.py 3``` for generating 10 subsets of data with 3 strains in each sample, followed by the execution of Demixer on each subset. 
The intermediate files will be saved in Demixer_cpp folder. Similarly for 4 strain data, run ```python simulate_dataset.py 4```.

Run ```python simulate_dataset_ablation.py 3``` for performing ablation analysis by generating 10 subsets of data with 3 strains in each sample.
Similarly for 4 strain data, run ```python simulate_dataset_ablation.py 4```.

### Simulation 2
For generating the reference strains of TBprofiler, run ```python generate_references.py```. Similarly run ```python generate_quanttb_references.py``` for generating the reference strains of QuantTB.

For generating the 8 different subsets consisting of 200 samples, run ```./exec_order.sh```. 

Run ```./deep_training_dataset.sh``` for generating artificial reads for training DeepMicrobes and DNABERT-S and then execute the below shell script for extracting the reads with mutations.

```
for i in *.sra_1.fastq; do
id=$(echo $i| cut -d "." -f 1);
bcftools view -v snps ${id}.vcf.gz |  bcftools filter -i 'QUAL>=20' > snps_only.vcf ;
convert2bed  -i vcf < snps_only.vcf > data.bed;
bedtools intersect -a ${id}.sorted_dedup_reads.bam -b data.bed -wa > output_reads.bam;
samtools view -h output_reads.bam > output_reads.sam;
samtools fastq output_reads.sam > output.fastq;
seqtk seq -A  output.fastq > ${id}.fasta;
rm snps_only.vcf; rm data.bed; rm output_reads.bam; rm output_reads.sam; rm output.fastq;
done

```

### Simulation 3

```for i in {0..11}; do python artsim_data_quant_split.py 7 5 5 $i; done``` - Run 10 iterations of Demixer on ARTmix data with the number of known strains set to 5. (few strains known + few unknown mode)

```for i in {0..11}; do python artsim_data_quant_split.py 7 5 7 $i; done``` - Run 10 iterations of Demixer on ARTmix data with the number of known strains set to 7. (all strains known mode)

```python artsim_data_quant_split.py 2 0 0 1``` - Run Demixer (iteration number 1) on 2-strain data with the number of known strains set to 0. (used for comparing vanilla LDA with Demixer)

Here, the first, second and third commandline arguments correspond to the number of strains, suffix of output_folder and the number of known strains respectively.

Note that the commands have to be exexuted from within the respective directories.
