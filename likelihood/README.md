# Likelihood estimation of allelic preference

1) Preprocessing of SNPs: `preprocessing_snps.sh`

   Script downloads ASB data from [AlleleDB](http://alleledb.gersteinlab.org/download)

   Files: ASB.auto.v2.1.aug16.txt.tgz   accB.auto.v2.1.aug16.txt.tgz

2) Scoring of binding affinity: `pyprobound_scoring.sh`

   Dependency: [PyProBound](https://github.com/BussemakerLab/PyProBound)

   Usage: `python ppb_score.py [model_fit] [input_seq_file] [output_ddG_file]`

3) Bootstrapping: `bootstrapping.sh` `bootstrapping_CTCF.sh`

   Usage: `Rscript bootstrapping.R [TF] [ref_ddG_file] [alt_ddG_file][boot_output_file]`

   Log-Likelihood model: `likelihood_functions.R`

4) Generating figures: `figures.sh`
