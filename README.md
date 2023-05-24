#     Mapping the semi-nested community structure of 3D chromosome contact networks

The goal of this project is to analyze the 3D structure of chromosomes using Hi-C data and network analysis techniques. Specifically, we aim to map out the chromosome's actual folding hierarchy by treating the measured DNA-DNA interactions as a weighted network and extracting 3D communities using the generalized Louvain algorithm.

In this project, we analyze Hi-C data from the human GM12878 B-lymphoblastoid cell line. We obtained published data from Gene Expression Omnibus (accession number is `GSE63525`). We normalized the Hi-C data using Knight-Ruiz Matrix Balancing (KR) algorithm implemented in `gcMapExplorer`.  We present KR-normalized Hi-C data for chromosomes 3, 5, 10, and 22 in the folder `data/Hi-C_data_KR_norm_chr3-5-10-22`.

The Hi-C data was treated as a network. Using the generalized Louvain algorithm, we tune its resolution parameter to scan through the community size spectrum, from A/B compartments to topologically associated domains (TADs), and then construct a tree connecting these communities. The community partitions, as well as the partitions for the irreducible domains, are published in folder `data/3Dcommunities_GenLouvain_output_across_gamma`.

We encourage users to contact the developers should problems arise.

# METHODS
In this repository, we provide code for community detection approach and other methods we describe in Methods Section of the manuscript.

More specifically, we provide an implementation of the GenLouvain algorithm, a MATLAB script, for detecting communities in Hi-C networks.
With the provided code, users can apply the algorithm to their own KR-normalized Hi-C datasets to identify communities of various size scales.

We also provide Python script to calculate communities' nestedness based on formulae given in [A new measure of ecological network structure based on node overlap and segregation](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12395).
Nestedness is a property of networks that quantifies the degree to which nodes with smaller degrees are "nested" within nodes with larger degrees. By calculating nestedness, users can gain insights into the underlying structure of communities and the relationships between them.

To characterize communities, we calculate HMM's Enrichment Folds for each community. Here we provide the script to do it. Enrichment folds help users understand the significance of specific HMM states within a given sample (community). The script includes a hypergeometric test to determine the probability of observing a certain number of peaks of an HMM state in a community and calculates the states' enrichment folds.

## GenLouvain community detection algorithm
_Script location:_ ```methods/GenLouvain/```<br>
_Script name:_ ```genLouvain_community_detection.m```<br>

_Description:_ This MATLAB script detects communities in chromosomes' Hi-C data using the GenLouvain community partition method.
For details of implementation, see Manuscript's Section II.B "The GenLouvain algorithm -- detecting 3D communities in Hi-C data"

_I. How to use this script?_<br>
To use this code, first specify which Hi-C data (chromosome name) will be partitioned into communities.
```matlab
chr1 = 10; %choose the first chromosome
chr2 = 10; %choose the second chromosome (for interchromsomal Hi-C) or repeat the first one (for intrachromosomal Hi-C)
```

Next, ensure that the Hi-C data and chromosome length information are available in the specified paths:
```matlab
data_path = ['/MATLAB Drive/mapping2023bernenko/HiC_data/chr' num2str(chr1) '_dna_00per.mat']
info_path = ['/MATLAB Drive/mapping2023bernenko/HiC_data/chr_' num2str(chr1) '.data.info']
```
For example, data for chromosome 10 are in files ```chr10_dna_00per.mat``` and ```chr_10.data.info```, and both files are located in folder ```methods/GenLouvain/HiC_data```

Afterwards, set the parameters:<br>
```alpha```: The exponent of the decay of contact frequencies between node i and node j. Default value is 1<br>
```gamma```: The resolution parameter, which can be set to a single value (for example 0.6) or a range of values (example is below).
```matlab
alpha = 1; % FIX one parameter

for gamma = 0.6% [0.4, 0.5, 0.6, 0.7, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9] ```
```
Last, run the algorithm. It will save the data into folder ```methods/GenLouvain/gL_partitions```

_II. Output_<br>
The script saves the community assignments for nodes in a .mat file in the gL_partitions directory, with the naming format chr{chr1}_gamma{gamma_n}_alpha{alpha_n}.mat.<br>

The output file contains:<br>
**ROWS**<br>
representing nodes;<br>
**COLUMNS**<br>
_column#1_ contains nodes' community assignment,<br>
_column#2_ is a chromosome name,<br>
_column#3_ is a modularity of the partition (not normalized).

## Nestedness

_Script location:_ ```methods/nestedness```<br>
_Script name:_ ```nestedness_Nij_calculation```

_Description:_ This Python script calculates the nestedness coefficient between two communities in a network based on their shared nodes. It also computes the expected number of shared nodes, probabilities, and various coefficients.
For details, see Section II.C "Network nestedness" of the manuscript.<br>
The code that implements formulae 3--10 is exported from the module ```methods/nestedness/utils/nestedness.py``` 

_Input Parameters_<br>
```nodes_shared:``` Number of shared nodes between community i and community j.<br>
```k_i:``` Degree of community i.<br>
```k_j:``` Degree of community j.<br>
```n_tot_num_nodes:``` Total number of domains in the network.

_Output_<br>
The script prints the following information:<br>
1. Expected number of shared nodes.<br>
2. W coefficient.<br>
3. Omega coefficient.<br>
4. Nestedness coefficient, N.<br>
5. Probabilities (less, exact, greater).

## Hypergeometric test and HMM's enrichment folds
_Script location:_```methods/hypergeometric_test```<br>
_Script name:_```hypergeometric_test_calculation```<br>

_Description:_ This Python script calculates the hypergeometric test for a given population, sample, and number of successful items. It computes the cumulative probability, confidence interval, mean, and enrichment folds.

_Input Parameters_<br>
```N:``` The number of items in the population.<br>
```k:``` The number of items in the population that are classified as successes.<br>
```n:``` The number of items in the sample.<br>
```x:``` The number of items in the sample that are classified as successes.

_Output_<br>
The script prints the following information:<br>
1. Cumulative probability for drawing the given number of successful items.
2. Confidence interval for the given alpha level.
3. Hypergeometric mean.
4. Enrichment folds for the given case.

For details about method implementation, please, see Section II.D "Chromatin states and folds of enrichment"
# DATA, how to read the files
## Hi-C data
Code snippet in Python to read Hi-C data from file `./data/Hi-C_data_KR_norm_chr3-5-10-22/chr3_dna_00per.mat`

```python
import scipy.io as sio
chrname = "chr3"
file_path = "./specify/your/path/"
adj = sio.loadmat(file_path + chrname +"_dna_00per.mat")
adj = adj['data']
```

## 3D communities of the Hi-C network
Code snippet in Python to read from file `./data/3Dcommunities_GenLouvain_output_across_gamma/chr3_gL_output.csv` that stores the community partitions for the range of gamma values (resolution parameter) and corresponding irreducible domains.

```python
import pandas as pd
chrname = "chr3"
file_path = "./specify/your/path/"
genlouvain_df = pd.read_csv(file_path+chrname+"_gL_output.csv", index_col=0)
```

Description of columns in a file `chr3_gL_output.csv`: 

| Column name      | Description |
| ----------- | ----------- |
| Index | ID of a Hi-C bin, given in a consecutive order from 0 to last |
| Bin_ID | the same as index |
| Bin_size(100kb) | the length of a DNA segment length. Each Hi-C bin has length of 1 (= 100 000 base pairs) |
| 0.9, …, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3 | resolution parameter gamma for the GenLouvain algorithm. The column entries represent the community IDs for each bin_ID. If not given (NaN), then this bin was not assigned to a community |
| Dom_ID | Irreducible domains are numbered consecutively from 0 to the last domain. A domain is defined as a group of linearly adjacent Hi-C bins that are members of the same community for each tested gamma (0.9, …, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3).|
