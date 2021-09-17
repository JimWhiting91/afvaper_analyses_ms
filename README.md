# AF-vapeR: Simulations and Empirical Demonstrations
This repository houses all scripts necessary to perform SLiM simulations of Full Parallel, Multi-Parallel, Divergent, and Variable Recombination simulations outlined in the manuscript:

*AF-vapeR: A multivariate genome scan for detecting parallel evolution using allele frequency change vectors*
doi:XXX

### Author Contact Information
james.whiting@ucalgary.ca

### Usage and license information
If you use or are inspired by code in this repository please cite the following work or contact me about how to cite.

---
![Simulated Eigenvalue Peaks](./simmed_peaks.png?raw=true "Simulated Eigenvalue Peaks")

## SLiM Simulations
Scripts to perform and process slim simulations are housed in `/slim`.

The following SLiM scripts are used to simulate parallel evolution under different modes:
  * Full Parallel - `/slim/SGV_sweeps_all_parallel_trees_fullrecap.slim`
  * Multi Parallel - `/slim/SGV_sweeps_two_parallel_lineages_trees_fullrecap.slim`
  * Divergent - `/slim/SGV_sweeps_four_divergent_lineages_fullrecap.slim`
  * Variable Recombination - `/slim/SGV_sweeps_all_parallel_trees_full_recap_variable_recombination`

These simulations are run by using the `slim/afvaper_slim_master_multi_parameters.sh` script, which is a bash array script designed for SLURM systems that parses parameters to slim simulations and processes outputs for Fst and AF-vapeR analyses.

The python script `slim/mutate_parallel_sweep_trees_cmd.py` is used for recapitation and mutation of SLiM treesequences to produce VCFs for AF-vapeR and Fst.

The R script `afvape_slim_vcf.R` performs AF-vapeR over slim VCFs. Both this step and the conversion of treesequences to VCF are handled by the master bash script.

## Empirical Datasets
Scripts for analysis of empirical datasets are stored in `R/`. This directory also handles some general plotting scripts.

The processing of VCFs to AF-vapeR outputs is done by the three `R/afvaper_test_*.R` scripts.

The subsequent analysis and visualisation of these outputs is done by the four `R/*bring_windows_together*.R` scripts.

Calculation and analysis of FPR and FNR of simulated data is handled by `R/fpr_fnr_analysis_of_sim_results.R`

General plotting scripts:
  * `R/plot_peaks_under_different_recombination.R` - takes simulated outputs from Variable Recombination simulations and visualises parallelism peaks.
  * `R/3d_scatter.R` - produces 3D scatter seen in Figure 1.
  * `R/chr15_vector_lines.R` - produces Figure S10 based on HP-LP guppy tests.
  * `R/visualise_sim_eigenvalue_examples.R` - produces Figure 2 based on sim outputs.
