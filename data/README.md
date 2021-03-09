# Title of the paper

This is an authors' implementation of "Title of the paper" in Python.

Authors: Meijun Gao, Kevin J. Liu

Paper: https://

Accepted by 

## Introduction

This work presents a new phylogenetic hidden Markov model-based method to analyze the GC-biased gene conversion and recombination
hotspots in eukaryotic genomes.

The new HMM model allows for modeling substitution parameters variation, mutation variation, and recombination rate variation,
and can be trained to detect the recombination hotspots and estimate substitution parameters.

Detailed information about the new phylogenetic HMM model is provided in https://......

## Simulation and empirical data 
This directory includes all simulation data in our paper.
+ 4 taxa 

    Sequence length=5k; background recombination rate=5; background mutation rate=1.
    - one recombination hotspot, recombinaiton rate multiple=10, mutation multiple=10; `data\simulation_data\h1_r10m10`
    - two recombination hotspots, recombinaiton rate multiple=10, mutation multiple=10; `data\simulation_data\h2_r10m10`
    - none recombination hotspot, negative control; use background recombination rate and mutation rate. `data\simulation_data\h0`
+ 5 taxa

    Sequence length=2k; background recombination rate=2; background mutation rate=1.
    - one recombination hotspot, recombinaiton rate multiple=10, mutation multiple=10; `data\simulation_data\h1_r10m10`
    - two recombination hotspots, recombinaiton rate multiple=10, mutation multiple=10; `data\simulation_data\h2_r10m10`
    - none recombination hotspot, negative control; use background recombination rate and mutation rate.  `data\simulation_data\h0`

The usage of `generate_data.py`:
```
python generate_hotspot.py --numtaxa=4 --modelname='HKY' --hotspot_num=1 
    --rr_multiple=10 --mr_multiple=10 --path=''
Options:
-h show the basic usage for simulation data generation
--numtaxa 4 or 5;   the length of the sequences are 5000, 2000 for 4 and 5 taxa, respectively
--modelname HKY or GTR; default HKY
--hotspot_num 0 or 1 or 2
    when hotspot_num=0, --rr_multiple --mr_multiple do not need  
--hot_rho recombination rate (rr) multiple based on the rr in background
--hot_dScale mutation rate (mr) multiple based on the mr in background
--path PATH, where to put output files
```
When you run the code, you need have the same directory structure, include `ms`, `msHOT`, `seq-gen`, `myFun.py`, 
and the directory `\tree_files` in the `simulation_src` directory. 

The previous two files can be downloaded in Hudson Lab Home Page[^2], `seq-gen` can be downloaded in [^3],
`myFun.py` is accessorial file. As for `\tree_files`, it includes all unrooted phylogenetic trees for two different
number of taxa (`\tree_files\taxa4_dif_topo.txt`, `\tree_files\taxa4_dif_topo.txt`) and random model trees to 
simulate the phylogenetic histories for simulation data (`\tree_files\taxa4_modeltree.txt`, `\tree_files\taxa5_modeltree.txt`).

[^2]: http://home.uchicago.edu/~rhudson1/
[^3]: https://snoweye.github.io/phyclust/document/Seq-Gen.v.1.3.2/Seq-Gen.Manual.html

### Empirical data
+ Simulation 

    Sequence length=5k; background recombination rate=5; background mutation rate=1.
    - one recombination hotspot, recombinaiton rate multiple=10, mutation multiple=10; `data\simulation_data\h1_r10m10`
    - two recombination hotspots, recombinaiton rate multiple=10, mutation multiple=10; `data\simulation_data\h2_r10m10`
    - none recombination hotspot, negative control; use background recombination rate and mutation rate. `data\simulation_data\h0`
+ 5 taxa

    Sequence length

