#!/usr/bin/env bash

#SNPs
gatk SelectVariants --variant results/variants/trio.vcf --select-type-to-include SNP --output results/variants/trio.SNP.vcf

gatk VariantFiltration --variant trio.SNP.vcf \
--filter-expression "QD < 2.0"              --filter-name "QD2" \
--filter-expression "QUAL < 30.0"           --filter-name "QUAL30" \
--filter-expression "SOR > 3.0"             --filter-name "SOR3" \
--filter-expression "FS > 60.0"             --filter-name "FS60" \
--filter-expression "MQ < 40.0"             --filter-name "MQ40" \
--filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
--output trio.SNP.filtered.vcf


#INDELs
gatk SelectVariants --variant results/variants/trio.vcf --select-type-to-include INDEL --output results/variants/trio.INDEL.vcf

gatk VariantFiltration --variant trio.INDEL.vcf \
--filter-expression "QD < 2.0"                  --filter-name "QD2" \
--filter-expression "QUAL < 30.0"               --filter-name "QUAL30" \
--filter-expression "FS > 200.0"                --filter-name "FS200" \
--filter-expression "ReadPosRankSum < -20.0"    --filter-name "ReadPosRankSum-20" \
--output trio.INDEL.filtered.vcf

#Merging filtered files
gatk MergeVcfs --INPUT trio.SNP.filtered.vcf --INPUT trio.INDEL.filtered.vcf --OUTPUT trio.filtered.vcf