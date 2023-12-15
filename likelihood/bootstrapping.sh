#!/usr/bin/env bash

mkdir data/bootstrap

### Control models
TFs=( CTCF EBF1 SPI1 )
for tf in "${TFs[@]}";
do
  Rscript bootstrapping_ctrl.R ${tf} data/bootstrap/${tf}_control_boot.rds
done

### MotifCentral models
IDs=( 12715 13774 13165 )
for i in "${!TFs[@]}";
do
	Rscript bootstrapping.R ${TFs[i]} \
    data/scores/${TFs[i]}_fit${IDs[i]}_k30_ref.tsv \
    data/scores/${TFs[i]}_fit${IDs[i]}_k30_alt.tsv \
    data/bootstrap/${TFs[i]}_fit${IDs[i]}_k30_boot.rds
done

### External resources
resourceIDs=(1291 2116)
for resourceID in "${resourceIDs[@]}";
do
	Rscript bootstrapping.R CTCF \
		data/scores/CTCF_external${resourceID}_k30_ref.tsv \
	  data/scores/CTCF_external${resourceID}_k30_alt.tsv \
    data/bootstrap/CTCF_external${resourceID}_k30_boot.rds
done

### PyProBound models
models=(
	CUTnTag_center   CUTnTag_end0
	ChIP-exo5_center ChIP-exo5_end1
	ChIP-seq_center  ChIP-seq_end5
)
for model in "${models[@]}";
do
	Rscript bootstrapping.R CTCF \
    data/scores/CTCF_${model}_k30_ref.tsv \
    data/scores/CTCF_${model}_k30_alt.tsv \
    data/bootstrap/CTCF_${model}_k30_boot.rds
done
