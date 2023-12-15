#!/usr/bin/env bash

mkdir data/scores
python import.py

### MotifCentral models
TFs=( CTCF EBF1 SPI1 )
motifCentralIDs=( 12715 13774 13165 )
for i in "${!TFs[@]}";
do
	python ppb_score.py data/bindingModels/motifcentral_fit_${motifCentralIDs[i]}.pt \
		data/seq/${TFs[i]}_k30_ref.txt \
		data/scores/${TFs[i]}_fit${motifCentralIDs[i]}_k30_ref.tsv
	python ppb_score.py data/bindingModels/motifcentral_fit_${motifCentralIDs[i]}.pt \
		data/seq/${TFs[i]}_k30_alt.txt \
		data/scores/${TFs[i]}_fit${motifCentralIDs[i]}_k30_alt.tsv
done

### External resources
resources=(HOCOMOCO_CTCF_HUMAN.H11MO.0.A JASPAR_CTCF_MA0139)
resourceIDs=(1291 2116)
for i in "${!resources[@]}";
do
	python ppb_score.py data/bindingModels/${resources[i]}.pt \
		data/seq/CTCF_k30_ref.txt \
		data/scores/CTCF_external${resourceIDs[i]}_k30_ref.tsv
	python ppb_score.py data/bindingModels/${resources[i]}.pt \
		data/seq/CTCF_k30_alt.txt \
		data/scores/CTCF_external${resourceIDs[i]}_k30_alt.tsv
done

### PyProBound models
models=(
	CUTnTag_center   CUTnTag_end0
	ChIP-exo5_center ChIP-exo5_end1
	ChIP-seq_center  ChIP-seq_end5
)
for model in "${models[@]}";
do
	python ppb_score.py ../in-vivo/fits/CTCF_${model}-ctcf_psam.pt \
		data/seq/CTCF_k30_ref.txt data/scores/CTCF_${model}_k30_ref.tsv
	python ppb_score.py ../in-vivo/fits/CTCF_${model}-ctcf_psam.pt \
		data/seq/CTCF_k30_alt.txt data/scores/CTCF_${model}_k30_alt.tsv
done
