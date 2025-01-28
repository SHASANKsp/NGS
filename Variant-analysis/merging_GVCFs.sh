#!/usr/bin/env bash


#combining GVCFs for different sample
gatk GenomicsDBImport --variant results/variants/mother.HC.g.vcf --variant results/variants/father.HC.g.vcf \
--variant results/variants/son.HC.g.vcf --intervals chr20:10018000-10220000 --genomicsdb-workspace-path results/genomicsdb

gatk GenotypeGVCFs --reference data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
--variant gendb://results/genomicsdb --intervals chr20:10018000-10220000 --output results/variants/trio.vcf


