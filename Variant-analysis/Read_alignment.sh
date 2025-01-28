#!/usr/bin/env bash

#Indexing reference genome with bwa
bwa index Homo_sapiens.GRCh38.dna.chromosome.20.fa

#alignment
bwa mem data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa data/fastq/mother_R1.fastq.gz data/fastq/mother_R2.fastq.gz > results/alignments/mother.sam

#alignment stats
samtools flagstat mother.sam > mother.sam.flagstat

#alignment sorting based on coordinate
samtools sort -o alignments/mother.sorted.sam alignments/mother.sam
 
#compressing sam to bam
samtools view -bh alignments/mother.sorted.sam > alignments/mother.bam

#Adding readgroups
gatk AddOrReplaceReadGroups \
--INPUT alignments/mother.bam \
--OUTPUT alignments/mother.rg.bam \
--RGLB lib1 \
--RGPU H0164.2.ALXX140820 \
--RGPL ILLUMINA \
--RGSM mother \
--RGID H0164.2.mother

#Marking of duplicates
gatk MarkDuplicates \
--INPUT alignments/mother.rg.bam \
--OUTPUT alignments/mother.rg.md.bam \
--METRICS_FILE alignments/marked_dup_metrics_mother.txt 

#alignment stats
samtools flagstat mother.rg.md.bam > mother.rg.md.bam.flagstat

#Indexing
samtools index mother.rg.md.bam

