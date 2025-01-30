#!/usr/bin/env bash


#SRA-TOOLKIT
#paired-end
./fasterq-dump SRR___ --split-files

#fastqc

#triming adapter sequences
java -jar trimmomatic-0.39.jar PE SRRforward.fastq SRRreverse.fastq outputforwardpaired.fastq.gz outputforwardunpaired.fastq.gz outputreversepaired.fastq.gz outputreverseoutunpaired.fastq.gz ILLUMINACLIP: path-to-adapter adapterUsed:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#Reference based alignment
#Indexing the refrence Genome
bowtie2-build MtbComplete_genome.fasta bowtie2
#alignment
bowtie2 -x path-to-index/bowtie2 -1 SRRfwdPaired.fastq -2  SRRreversePaired.fastq -S path-for-outputfile/OUTsam

#processing
samtools view -S -b OUTsam > OUT.bam
samtools sort OUT.bam -o OUTsorted.bam
samtools index OUTsorted.bam #Checkpoint-1

#Separating Illumina Based Strand Specific Rna-Seq Alignments 
DATA = "OUTsorted.bam"

# Forward strand.
samtools view -b -f 128 -F 16 $DATA > fwd1.bam
samtools index fwd1.bam
samtools view -b -f 80 $DATA > fwd2.bam
samtools index fwd2.bam

# Combine alignments that originate on the forward strand.
samtools merge -f fwd.bam fwd1.bam fwd2.bam
samtools index fwd.bam

# Reverse strand
samtools view -b -f 144 $DATA > rev1.bam
samtools index rev1.bam
samtools view -b -f 64 -F 16 $DATA > rev2.bam
samtools index rev2.bam

# Combine alignments that originate on the reverse strand.
samtools merge -f rev.bam rev1.bam rev2.bam
samtools index rev.bam


#Convert GFF to BED file
convert2bed --input=gff --output=bed  <mtb_Genome.gff>  mtb_completeGenome.bed

#Map Bam to Bed file
bedtools multicov -bams fwd.bam -bed  mtb_completeGenome.bed > fwd.bed
bedtools multicov -bams rev.bam -bed  mtb_completeGenome.bed > rev.bed

#Per Base Coverage Bed 
Forward: bedtools genomecov -ibam fwd.bam -d > SRR21675925perbaseFWD.bed 
Reverse: bedtools genomecov -ibam rev.bam -d > SRR21675925perbaseREV.bed