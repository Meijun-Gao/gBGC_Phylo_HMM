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

## Simulation data 
This directory includes all simulation data in our paper. Here list all simulation conditions:
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

In each condition directory, it includes 20 `REPLICATE` directories, a `new` directory and a `simple` directory. Each 
`REPLICATE` directory includes a `sequences.fasta` file and a `all-genetrees` file that contains all possible unrooted 
topology trees with given number of taxa.

The `new` directory includes the learning output files of the new Phylo-HMM based algorithm.

The `simple` directory includes the learning output files of the simple Phylo-HMM based algorithm.

## Empirical data
The `empirical\wholeseq\` includes the concatenated aligned gene sequences of 12 chromosomes. 

The `empirical\SNPseq\` includes the corresponding SNP sequences in `empirical\wholeseq\`.

The `result` includes the learning output files of the new and simple Phylo-HMM based algorithm in each subdirectory.

