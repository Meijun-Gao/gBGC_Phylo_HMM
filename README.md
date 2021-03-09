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

## Requirements

python == 3.7

Set up the conda env: conda env create -f hmm_env.yml

## Run the model

To get the posterior decoding results of all states and estimated parameters, run:
```
python run_HMM.py alignmentFile.fasta all-genetrees.txt 
Options:
    -h show the basic usage of the algorithm
    -t INT number of independent trials to run;  default value=10
    --sh_ac sh activate; default value=1
        1: background and hotspot use different substitution model parameters, Mb!=Mh
        0: background and hotspot use different substitution model parameters, Mb=Mh
    --modelname substitutiom model name; HKY or GTR;  default value='HKY'
    --prefix PATH, where to put output files
    --SNP if use the SNP sequences of the input sequences; default value=0
        0 use all site;  1 only use SNP
```
The `\RAXML` directory is the complied RAxML files, which can be downloaded from the [^1]

[^1]: https://cme.h-its.org/exelixis/web/software/raxml/
## Output file explanations

The output file of the algorithm includes:
+ Optimize time: xxx hours
+ Begin Prob: the log-likelihood value at the beginning of the learning process
+ End Prob: the log-likelihood value after all learning process
+ Substitution parameters back: the substitution model parameters for background regions in genome sequences.
    - format: (f_A, f_C, f_G, f_T, ts/tv), the first four parameters represent the base frequencies of A, C, G, and T, ts/tv represents
the ratio of transition and transversion substitutions.
+ Substitution parameters hot: the substitution model parameters for hotspot regions in genome sequences.
    - format: the same as before
+ Transition parameters: the parameters for the transition matrix in the paper, (sb, sh, $\gamma$)
+ hscale parameters: the scaling coefficient for hotspot genealogy tree branch lengths, $\beta$
+ Trees: list all gene trees corresponding to each state
+ logarithmic posterior probability of each state in each site of the input sequences
    - format: Position, Posteriors1, Posteriors2, Posteriors3, ..., Posteriors(2K)


## Simulation and empirical data explanations
We use msHOT and seq-gen to generate simulation data. You can use `generate_data.py` to generate simulation data.
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
The empirical data include 4 rice genomes' gene alignments, Oryza Sativa Japonica, Oryza Sativa Indica, Oryza Rufipogon,
and Oryza Nivara. We obtain their orthologous genes from EnsemblPlants[^4], align these genes by MAFFT[^5], and then concatenate
them in chromosome level according to their gene annotations. We obtained 12 sequences files (`\empirical\wholeseq\`), 
each of them includes four aligned sequence corresponding to the four rice. In order to speed up the running time, we 
use the SNP sequences (`\empirical\SNPseq\`) as the input, although the code can transfer the whole sequences into SNP sequences.   

[^4]: http://plants.ensembl.org/index.html
[^5]: https://mafft.cbrc.jp/alignment/software/
## Citation

If using this code, please cite our paper.
```

```# gBGC_Phylo_HMM
# gBGC_Phylo_HMM