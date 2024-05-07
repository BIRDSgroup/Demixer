# Demixer: A multi-sample analysis tool to identify mixed infection

This is the official repository of the Manuscript "Demixer: A probabilistic generative model to delineate different strains of a microbial species in a mixed infection sample" by Brintha V P, and Manikandan Narayanan. Please note that the [Figures](https://github.com/BIRDSgroup/Demixer/tree/main/Figures) folder has the codes for reproducing the figures in the manuscript. 

## License preamble 

Copyright 2024 BIRDS Group, IIT Madras

Demixer is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Demixer is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with Demixer.  If not, see <https://www.gnu.org/licenses/>.

## Installation

See ```requirements.txt``` file for the list of dependencies. The main program Demixer is implemented using C++, and written in compliance with the C++11 standard. Python scripts used for pre-processing of the data are also provided. 

Download all the required files using the below command.
```
git clone https://github.com/BIRDSgroup/Demixer.git
```

## Usage

### Running Demixer on a set of samples using .vcf input
Demixer takes a .vcf file as input and processes it to generate the Sample-SNP matrix. For preprocessing, the below command has to be run. It takes 4 input parameters. 
- input_.vcf_file (multisample .vcf file)
- output_folder (the folder name in which the intermediate files are to be saved)
- db_name (the reference database name: tbprof, quanttb or covid)
- AD_COV  (.vcf type: 2 [for calldata/AD] or 3[for calldata/COV]).

The output files resulting from preprocessing will be saved in the finaloutput/output_folder. The sample-SNP matrix will be saved in Docs_int_onlySnpgap.dat file. By default, the new folder output_folder will be created within the finaloutput folder.
```
python main_preprocess.py input_.vcf_file output_folder db_name AD_COV
```
Once the preprocessing is done, run the CGS algorithm on the sample-SNP matrix using the executable file Demixer generated from the main.cpp file. It takes one input, the output_folder containing the intermediate files. The model parameters $\theta$ of size $m$ x $k$ and $\phi$ of size $k$ x $n$ will be saved in n_m_z0.dat and n_z_t0.dat files respectively in the output_folder directory. The number of strains $k$ will be determined during preprocessing, $m$ represents the number of samples and $n$ represents the number of SNP-alleles. The proportions of the strains in the samples can be viewed in the finaloutput/output_folder/n_m_z0.txt file.
```
./Demixer finaloutput/output_folder
```
Run the command ```chmod +x Demixer```, if running the above command for the first time.

The folder finaloutput has the subfolders corresponding to each run of Demixer on different .vcf files. After the running of CGS algorithm, the inferred parameters need to be postprocessed by running the following command. The third parameter takes the value True/False and by default it is set to **False** (no merging of 2 or more *de novo* strains). The output files resulting from postprocessing will be saved in the finaloutput/output_folder. The inferred proportions, lineage names, confidence and mixed/non-mixed call for each sample are saved in the files prop.csv, lineage.csv, confidence_new.csv  and Mixpred.csv respectively. Also the mode ($o_m$) and frequency of mode ($f_m$) for each sample determined using the SNP plot is saved in the file snp_plot_conf.csv respectively.

```
python postprocessing.py finaloutput/output_folder/
python postprocessing.py finaloutput/output_folder/ True (allows merging of *de novo* strains)    
```

The postprocessing results reported in the manuscript are obtained by running the above command using the default value for the third parameter.

### Testing new samples in .vcf file using CRyPTIC-trained parameters
To test new samples, run the below commands to generate the sample-SNP matrix, run CGS algorithm and postprocessing. Requires 3 input parameters: input_.vcf_file, db_name and AD_COV. The intermediate and output files will be saved in finaloutput/test directory. A sample .vcf file for testing and the trained $\phi$ parameter can be downloaded from this [link](https://drive.google.com/drive/folders/1-zzEhnMofpfUvxH17KaJ23qp_SJzu9R4?usp=drive_link). Move the n_z_t0.dat file ($\phi$ parameter) to *scripts/finaloutput/Cryptic_rerun* folder. For running this file set db_name to **tbprof** and AD_COV to **2**.
```
python test.py input_.vcf_file test db_name AD_COV
./test_Demixer;
python test_strain_ids.py finaloutput/test/
```
Run ```chmod +x test_Demixer```, if running the above command for the first time. Also, if there are any issues while running test_Demixer, please compile the test.cpp file using the commands in the next section. Make sure that the above commands are ran from within the scripts folder.

### Generating executable code from .cpp file
For generating the Demixer (training mode) executable file, the following command is run:

```
g++ -g -O4  main.cpp globals.cpp -std=c++11 -Wall -ggdb -pg -g3 -fopenmp -lgsl -lgslcblas  -funroll-loops -pg -Ofast -o Demixer
```

For generating the test_Demixer (testing mode) executable file, the following command is run:

```
g++ -g -O4  test.cpp globals.cpp -std=c++11 -Wall -ggdb -pg -g3 -fopenmp -lgsl -lgslcblas  -funroll-loops -pg -Ofast -o test_Demixer
```
Note that all the relavant codes are in the scripts folder and the above commands should be executed from within the scripts directory.

## Folder Description

The [Figures](https://github.com/BIRDSgroup/Demixer/tree/main/Figures) folder has the codes for reproducing the figures in the manuscript. 

The [scripts](https://github.com/BIRDSgroup/Demixer/tree/main/scripts) folder has the codes for running Demixer on a set of samples from .vcf file.

The [simulation](https://github.com/BIRDSgroup/Demixer/tree/main/simulation) folder has the scripts and other files required for generating simulated data and running Demixer on those datasets.

The [scripts/db](https://github.com/BIRDSgroup/Demixer/tree/main/scripts/db) folder has the reference database files related to TBprofiler and Quanttb.

The [scripts/finaloutput](https://github.com/BIRDSgroup/Demixer/tree/main/scripts/finaloutput) folder has the sub-folders containing the intermediate files related to each run. 

