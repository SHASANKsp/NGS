#!/usr/bin/env bash


cellranger count --id=ETV6-RUNX1_1 --sample=ETV6-RUNX1_1 --transcriptome=cellranger_index  --fastqs=course_data/reads --localcores=4 --create-bam=true