# MorphoModels
 <p align="center">
  <img src="https://github.com/sergeitarasov/MorphoModels/blob/main/vignettes/Fig_icon.png" width="700" title="hover text">
</p>  

# Content

This repository contains [vignette](https://github.com/sergeitarasov/MorphoModels/wiki/Constructing-rate-matrices-for-ED-models) and **scripts** for setting up and running Embedded Dependency (ED) models for phylogenetic inference with inapplicable morphological data [(Tarasov 2022)](https://www.biorxiv.org/content/10.1101/2021.04.26.441495v3).

* [Here is a vignette](https://github.com/sergeitarasov/MorphoModels/wiki/Constructing-rate-matrices-for-ED-models) for setting up ED models and rate matrices to run in RevaBayes or other R packages.
* Below are the **scripts** for reproducing analyses from [Tarasov (2022)](https://www.biorxiv.org/content/10.1101/2021.04.26.441495v3).
 
## The scripts used in [Tarasov (2022)](https://www.biorxiv.org/content/10.1101/2021.04.26.441495v3)

You can download the entire repository by clicking on the green `Code` button and then selecting `Download ZIP` in the upper right corner. It may also be downloaded by using this code:

``` r
$ git clone https://github.com/sergeitarasov/MorphoModels
```

- **1_ED_model_performance** The scripts for "1. ED models and SMMs: fixed data".
- - RevBayes scripts for running tree inference with the ED models and calculating marginal likelihood under fixed data.

- **2_EDvsHMM** The scripts for "1. ED models and SMMs: fixed data".
- - RevBayes scripts for tree inference (and marginal likelihood calculation) with expanded ED models and SMM  under fixed data. 

- **3_ED_TailColor_problem** The scripts for "3. ED and Tail Color Problem".
- - RevBayes scripts to test the resolution of LTC in the tail color problem under fixed data. 

- **4a_EDvsHMM_topology** The scripts for "2. ED models vs. SMMs: simulated data".
- - R scripts (in *RB_run_1/4a_R*) to simulate trees and characters under the ED6 model.
- - RevBayes scripts to perform tree inferences under the ED6, ED3, and two SMMs.
- - R scripts (in *RB_run_1/4a_R*) to calculate boostrapped p-values.


# References

Tarasov, Sergei. 2022. “New Phylogenetic Markov Models for Inapplicable
Morphological Characters.” *bioRxiv*, 2021–04.
<https://www.biorxiv.org/content/10.1101/2021.04.26.441495v3.abstract>.
