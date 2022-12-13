# MorphoModels
 
This repository contains scripts for reproducing analyses from the study "Tarasov S., 2022. New phylogenetic Markov models for inapplicable morphological characters" available at [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.04.26.441495v2).
 
- **1_ED_model_performance** The scripts for "1. ED models and SMMs: fixed data".
- - RevBayes scripts for running tree inference with the ED models and calculating marginal likelihood under fixed data.

- **2_EDvsHMM** The scripts for "1. ED models and SMMs: fixed data".
- - RevBayes scripts for tree inference (and marginal likelihood calculation) with the expanded ED models and SMM  under fixed data. 

- **3_ED_TailColor_problem** The scripts for "3. ED and Tail Color Problem".
- - RevBayes scripts to test the resolution of LTC in the tail color problem under fixed data. 

- **4a_EDvsHMM_topology** The scripts for "2. ED models vs. SMMs: simulated data".
- - R scripts (in *RB_run_1/4a_R*) to simulate trees and characters under the ED6 model.
- - RevBayes scripts to perform tree inferences under the ED6, ED3, and two SMMs.
- - R scripts (in *RB_run_1/4a_R*) to calculate boostrapped p-values.

 <p align="left">
  <img src="https://github.com/sergeitarasov/MorphoModels/blob/main/vignettes/Fig_icon.png" width="100" title="hover text">
</p>  
