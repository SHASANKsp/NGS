#!/usr/bin/env bash

cd ~/workdir/alignment_output

samtools view -bh -F 4 SRR519926.sorted.bam > SRR519926.sorted.mapped.bam
samtools view -c SRR519926.sorted.mapped.bam
samtools view -bh -f 0x4 SRR519926.sorted.bam > SRR519926.sorted.unmapped.bam
samtools view -c SRR519926.sorted.unmapped.bam