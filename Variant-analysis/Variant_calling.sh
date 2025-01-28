#!/usr/bin/env bash

#Creating the indices for the vcfs
gatk IndexFeatureFile --input variants/1000g_gold_standard.indels.filtered.vcf
gatk IndexFeatureFile --input variants/GCF.38.filtered.renamed.vcf

#Creating the index for the reference genome
samtools faidx reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa
gatk CreateSequenceDictionary --REFERENCE reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa


#BQSR
gatk BaseRecalibrator --reference data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa --input results/alignments/mother.rg.md.bam \
--known-sites data/variants/GCF.38.filtered.renamed.vcf --known-sites data/variants/1000g_gold_standard.indels.filtered.vcf \
--output results/bqsr/mother.recal.table

gatk ApplyBQSR --input results/alignments/mother.rg.md.bam --bqsr-recal-file results/bqsr/mother.recal.table --output results/bqsr/mother.recal.bam

#variant calling
gatk HaplotypeCaller --reference data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
--input results/bqsr/mother.recal.bam --output results/variants/mother.HC.vcf --intervals chr20:10018000-10220000

gatk VariantsToTable --variant results/variants/mother.HC.vcf --fields CHROM -F POS -F TYPE -GF GT --output results/variants/mother.HC.table
