#!/usr/bin/env bash

mkdir data
mkdir data/accB
mkdir data/seq

wget -P data/ http://alleledb.gersteinlab.org/download/ASB.auto.v2.1.aug16.txt.tgz
wget -P data/ http://alleledb.gersteinlab.org/download/accB.auto.v2.1.aug16.txt.tgz

tar -xvzf data/ASB.auto.v2.1.aug16.txt.tgz -C data/
tar -xvzf data/accB.auto.v2.1.aug16.txt.tgz -C data/

Rscript preprocessing_snps.R