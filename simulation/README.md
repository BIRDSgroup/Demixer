## Simulation

This folder has the scripts used for generating the different synthetic data and running Demixer on them.

### Simulation 1
Run ```python simulate_dataset.py 3``` for generating 10 subsets of data with 3 strains in each sample, followed by the execution of Demixer on each subset. 
The intermediate files will be saved in Demixer_cpp folder. Similarly for 4 strain data, run ```python simulate_dataset.py 4```.

### Simulation 2
For generating the reference strains of TBprofiler, run ```python generate_references.py```. Similarly run ```python generate_quanttb_references.py``` for generating the reference strains of QuantTB.

For generating the 8 different subsets consisting of 200 samples, run ```./exec_order.sh```. The shell script simulation2run.sh has to be run for the execution of Demixer on the different subsets of data.

### Simulation 3

```python artsim_data_quant_split.py 7 5 5``` - Run Demixer on ARTmix_1 data with the number of known strains set to 5. (few strains known + few unknown mode)

```python artsim_data_quant_split.py 7 7 7``` - Run Demixer on ARTmix_1 data with the number of known strains set to 7. (all strains known mode)

```python artsim_data_quant_split.py 2 0 0``` - Run Demixer on ARTmix_2 data with the number of known strains set to 0. (used for comparing vanilla LDA with Demixer)

Here, the first, second and third commandline arguments correspond to the number of strains, suffix of output_folder and the number of known strains respectively.

