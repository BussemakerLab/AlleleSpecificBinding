#!/usr/bin/env bash

mkdir processing
wget -P processing https://downloads.sourceforge.net/project/bbmap/BBMap_39.01.tar.gz
tar -xvzf processing/BBMap_39.01.tar.gz -C processing/

mkdir processing/GRCh38_full_analysis_set
wget -P processing/GRCh38_full_analysis_set \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
gunzip processing/GRCh38_full_analysis_set/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz

wget -P processing https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.7/sratoolkit.3.0.7-mac64.tar.gz
tar -xvzf processing/sratoolkit.3.0.7-mac64.tar.gz -C processing/

fasterqdump=processing/sratoolkit.3.0.7-mac64/bin/fasterq-dump
bbmap=processing/bbmap
assembly=processing/GRCh38_full_analysis_set/GCA_000001405.15_GRCh38_full_analysis_set.fna

$bbmap/bbmap.sh ref=$assembly -Xmx40G

bam2bed() {
    samtools view --threads 6 -h -q 30 -F 1804 -f 2 -b "$1" \
        | samtools sort --threads 6 -m 20G -n - \
        | samtools fixmate - - \
        | bedtools bamtobed -bedpe -i stdin \
        | awk -v FS="\t" -v OFS="\t" '{split($1,c," "); print c[1], $2, $6}'
}

bam2bed_center() {
    samtools view --threads 6 -h -q 30 -F 1804 -f 2 -b "$1" \
        | samtools sort --threads 6 -m 20G -n - \
        | samtools fixmate - - \
        | bedtools bamtobed -bedpe -i stdin \
        | awk -v FS="\t" -v OFS="\t" '{split($1,c," "); a=int(($2+$6)/2+0.5); if (a>=100) print c[1], a-100, a+100}'
}

bed2fa() {
    bedtools getfasta -fi $assembly -bed "$1"
}

fa2tsv() {
    gunzip -c "$1" | grep -v \> | sort | uniq -c \
    | sed -e 's/ *//' -e 's/ /\t/' \
    | awk ' { t = $1; $1 = $2; $2 = t; print; } ' | tr ' ' \\t | tr a-z A-Z
}

fq2tsv() {
    gunzip -c "$1" | sed -n 2~4p | sort | uniq -c \
    | sed -e 's/ *//' -e 's/ /\t/' \
    | awk ' { t = $1; $1 = $2; $2 = t; print; } ' | tr ' ' \\t | tr a-z A-Z
}

sampletsv() {
    python3 sample.py "$1" "$2" "$3"
}

# CTCF ChIP-seq
wget -P data/CTCF_ChIP-seq \
    https://www.encodeproject.org/files/ENCFF873ZZU/@@download/ENCFF873ZZU.fastq.gz
wget -P data/CTCF_ChIP-seq \
    https://www.encodeproject.org/files/ENCFF611URR/@@download/ENCFF611URR.fastq.gz
$bbmap/bbmerge.sh in1=data/CTCF_ChIP-seq/ENCFF873ZZU.fastq.gz \
    in2=data/CTCF_ChIP-seq/ENCFF611URR.fastq.gz \
    out=data/CTCF_ChIP-seq/SRR5112225.fastq.gz
fq2tsv data/CTCF_ChIP-seq/SRR5112225.fastq.gz \
    | gzip > data/CTCF_ChIP-seq/SRR5112225.tsv.gz
wget -P data/CTCF_ChIP-seq \
    https://www.encodeproject.org/files/ENCFF631JSV/@@download/ENCFF631JSV.fastq.gz
wget -P data/CTCF_ChIP-seq \
    https://www.encodeproject.org/files/ENCFF715KYL/@@download/ENCFF715KYL.fastq.gz
$bbmap/bbmerge.sh in1=data/CTCF_ChIP-seq/ENCFF631JSV.fastq.gz \
    in2=data/CTCF_ChIP-seq/ENCFF715KYL.fastq.gz \
    out=data/CTCF_ChIP-seq/SRR5111853.fastq.gz
fq2tsv data/CTCF_ChIP-seq/SRR5111853.fastq.gz \
    | gzip > data/CTCF_ChIP-seq/SRR5111853.tsv.gz
sampletsv data/CTCF_ChIP-seq/SRR5112225.tsv.gz \
    data/CTCF_ChIP-seq/SRR5111853.tsv.gz \
    data/CTCF_ChIP-seq/Control-R1_CTCF-R1.1000000.tsv.gz

# CTCF ChIP-seq-centered
wget -P data/CTCF_ChIP-seq \
    https://www.encodeproject.org/files/ENCFF533HJT/@@download/ENCFF533HJT.bam
bam2bed_center data/CTCF_ChIP-seq/ENCFF533HJT.bam \
    > data/CTCF_ChIP-seq/SRR5112225-center.bed
bed2fa data/CTCF_ChIP-seq/SRR5112225-center.bed \
    | gzip > data/CTCF_ChIP-seq/SRR5112225-center.fa.gz
gzip data/CTCF_ChIP-seq/SRR5112225-center.bed
fa2tsv data/CTCF_ChIP-seq/SRR5112225-center.fa.gz \
    | gzip > data/CTCF_ChIP-seq/SRR5112225-center.tsv.gz
wget -P data/CTCF_ChIP-seq \
    https://www.encodeproject.org/files/ENCFF595PUL/@@download/ENCFF595PUL.bam
bam2bed_center data/CTCF_ChIP-seq/ENCFF595PUL.bam \
    > data/CTCF_ChIP-seq/SRR5111853-center.bed
bed2fa data/CTCF_ChIP-seq/SRR5111853-center.bed \
    | gzip > data/CTCF_ChIP-seq/SRR5111853-center.fa.gz
gzip data/CTCF_ChIP-seq/SRR5111853-center.bed
fa2tsv data/CTCF_ChIP-seq/SRR5111853-center.fa.gz \
    | gzip > data/CTCF_ChIP-seq/SRR5111853-center.tsv.gz
sampletsv data/CTCF_ChIP-seq/SRR5112225-center.tsv.gz \
    data/CTCF_ChIP-seq/SRR5111853-center.tsv.gz \
    data/CTCF_ChIP-seq/Control-R1_CTCF-R1.center.1000000.tsv.gz

# CTCF ChIP-Exo 5.0
$fasterqdump SRR6736398 -s -Z | gzip > data/CTCF_ChIP-exo5/SRR6736398.fastq.gz
$bbmap/bbmap.sh in=data/CTCF_ChIP-exo5/SRR6736398.fastq.gz \
    out=data/CTCF_ChIP-exo5/SRR6736398.bam -Xmx40G interleaved=True
bam2bed data/CTCF_ChIP-exo5/SRR6736398.bam \
    > data/CTCF_ChIP-exo5/SRR6736398.bed
bed2fa data/CTCF_ChIP-exo5/SRR6736398.bed \
    | gzip > data/CTCF_ChIP-exo5/SRR6736398.fa.gz
gzip data/CTCF_ChIP-exo5/SRR6736398.bed
fa2tsv data/CTCF_ChIP-exo5/SRR6736398.fa.gz \
    | gzip > data/CTCF_ChIP-exo5/SRR6736398.tsv.gz
$fasterqdump SRR6736390 -s -Z | gzip > data/CTCF_ChIP-exo5/SRR6736390.fastq.gz
$bbmap/bbmap.sh in=data/CTCF_ChIP-exo5/SRR6736390.fastq.gz \
    out=data/CTCF_ChIP-exo5/SRR6736390.bam -Xmx40G interleaved=True
bam2bed data/CTCF_ChIP-exo5/SRR6736390.bam \
    > data/CTCF_ChIP-exo5/SRR6736390.bed
bed2fa data/CTCF_ChIP-exo5/SRR6736390.bed \
    | gzip > data/CTCF_ChIP-exo5/SRR6736390.fa.gz
gzip data/CTCF_ChIP-exo5/SRR6736390.bed
fa2tsv data/CTCF_ChIP-exo5/SRR6736390.fa.gz \
    | gzip > data/CTCF_ChIP-exo5/SRR6736390.tsv.gz
sampletsv data/CTCF_ChIP-exo5/SRR6736398.tsv.gz \
    data/CTCF_ChIP-exo5/SRR6736390.tsv.gz \
    data/CTCF_ChIP-exo5/noAb_CTCF.1000000.tsv.gz

# CTCF ChIP-Exo 5.0 centered
bam2bed_center data/CTCF_ChIP-exo5/SRR6736398.bam \
    > data/CTCF_ChIP-exo5/SRR6736398-center.bed
bed2fa data/CTCF_ChIP-exo5/SRR6736398-center.bed \
    | gzip > data/CTCF_ChIP-exo5/SRR6736398-center.fa.gz
gzip data/CTCF_ChIP-exo5/SRR6736398-center.bed
fa2tsv data/CTCF_ChIP-exo5/SRR6736398-center.fa.gz \
    | gzip > data/CTCF_ChIP-exo5/SRR6736398-center.tsv.gz
bam2bed_center data/CTCF_ChIP-exo5/SRR6736390.bam \
    > data/CTCF_ChIP-exo5/SRR6736390-center.bed
bed2fa data/CTCF_ChIP-exo5/SRR6736390-center.bed \
    | gzip > data/CTCF_ChIP-exo5/SRR6736390-center.fa.gz
gzip data/CTCF_ChIP-exo5/SRR6736390-center.bed
fa2tsv data/CTCF_ChIP-exo5/SRR6736390-center.fa.gz \
    | gzip > data/CTCF_ChIP-exo5/SRR6736390-center.tsv.gz
sampletsv data/CTCF_ChIP-exo5/SRR6736398-center.tsv.gz \
    data/CTCF_ChIP-exo5/SRR6736390-center.tsv.gz \
    data/CTCF_ChIP-exo5/noAb_CTCF.center.1000000.tsv.gz

# CTCF CUT&Tag
$fasterqdump SRR8435051 -s -Z | gzip > data/CTCF_ChIP-CUTnTag/SRR8435051.fastq.gz
$bbmap/bbmap.sh in=data/CTCF_CUTnTag/SRR8435051.fastq.gz \
    out=data/CTCF_CUTnTag/SRR8435051.bam -Xmx40G interleaved=true
bam2bed data/CTCF_CUTnTag/SRR8435051.bam \
    > data/CTCF_CUTnTag/SRR8435051.bed
bed2fa data/CTCF_CUTnTag/SRR8435051.bed \
    | gzip > data/CTCF_CUTnTag/SRR8435051.fa.gz
gzip data/CTCF_CUTnTag/SRR8435051.bed
fa2tsv data/CTCF_CUTnTag/SRR8435051.fa.gz \
    | gzip > data/CTCF_CUTnTag/SRR8435051.tsv.gz
$fasterqdump SRR8754587 -s -Z | gzip > data/CTCF_ChIP-CUTnTag/SRR8754587.fastq.gz
$bbmap/bbmap.sh in=data/CTCF_CUTnTag/SRR8754587.fastq.gz \
    out=data/CTCF_CUTnTag/SRR8754587.bam -Xmx40G interleaved=true
bam2bed data/CTCF_CUTnTag/SRR8754587.bam \
    > data/CTCF_CUTnTag/SRR8754587.bed
bed2fa data/CTCF_CUTnTag/SRR8754587.bed \
    | gzip > data/CTCF_CUTnTag/SRR8754587.fa.gz
gzip data/CTCF_CUTnTag/SRR8754587.bed
fa2tsv data/CTCF_CUTnTag/SRR8754587.fa.gz \
    | gzip > data/CTCF_CUTnTag/SRR8754587.tsv.gz
sampletsv data/CTCF_CUTnTag/SRR8435051.tsv.gz \
    data/CTCF_CUTnTag/SRR8754587.tsv.gz \
    data/CTCF_CUTnTag/IgG_CTCF.1000000.tsv.gz

# CTCF CUT&Tag centered
bam2bed_center data/CTCF_CUTnTag/SRR8435051.bam \
    > data/CTCF_CUTnTag/SRR8435051-center.bed
bed2fa data/CTCF_CUTnTag/SRR8435051-center.bed \
    | gzip > data/CTCF_CUTnTag/SRR8435051-center.fa.gz
gzip data/CTCF_CUTnTag/SRR8435051-center.bed
fa2tsv data/CTCF_CUTnTag/SRR8435051-center.fa.gz \
    | gzip > data/CTCF_CUTnTag/SRR8435051-center.tsv.gz
bam2bed_center data/CTCF_CUTnTag/SRR8754587.bam \
    > data/CTCF_CUTnTag/SRR8754587-center.bed
bed2fa data/CTCF_CUTnTag/SRR8754587-center.bed \
    | gzip > data/CTCF_CUTnTag/SRR8754587-center.fa.gz
gzip data/CTCF_CUTnTag/SRR8754587-center.bed
fa2tsv data/CTCF_CUTnTag/SRR8754587-center.fa.gz \
    | gzip > data/CTCF_CUTnTag/SRR8754587-center.tsv.gz
sampletsv data/CTCF_CUTnTag/SRR8435051-center.tsv.gz \
    data/CTCF_CUTnTag/SRR8754587-center.tsv.gz \
    data/CTCF_CUTnTag/IgG_CTCF.center.1000000.tsv.gz
