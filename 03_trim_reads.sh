#!/usr/bin/env bash

TRIMMED_DIR=~/workdir/trimmed_data
READS_DIR=~/workdir/reads

mkdir -p $TRIMMED_DIR

cutadapt \
--adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
[QUALITY CUTOFF OPTION] \
[MINIMUM LENGTH OPTION] \
--output $TRIMMED_DIR/paired_trimmed_SRR519926_1.fastq \
--paired-output $TRIMMED_DIR/paired_trimmed_SRR519926_2.fastq \
$READS_DIR/SRR519926_1.fastq \
$READS_DIR/SRR519926_2.fastq
