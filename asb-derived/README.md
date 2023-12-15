# De novo motif discovery from allele-aware ChIP-seq counts

`asb.py` contains the implementation of the allelic preference likelihood model in
[PyProBound](https://github.com/BussemakerLab/PyProBound).

`run.sh` processes the ASB data, fits the ASB models, and generates the relevant figures.
It assumes that the likelihood scripts ([1](../likelihood/preprocessing_snps.sh),
[2](../likelihood/pyprobound_scoring.sh)) have already been run.
