#!/bin/bash

for file in /home/exx/Desktop/Jake/M32_GBS/trimmed_reads/*.fastq.gz; do

name="$(basename $file .fastq.gz)"
echo "working on $name"
bwa mem -R $(echo "@RG\tID:A00223_687_HM577DRXY_1\tPL:illumina\tPU:A00223_687_HM577DRXY_1\tLB:$name\tSM:$name") -t 85 /home/exx/Desktop/Jake/M32_GBS/genome/Bs_v2_soft.mask.fasta $file > /home/exx/Desktop/Jake/M32_GBS/aligned_reads/$name.sam
samtools view -S -b /home/exx/Desktop/Jake/M32_GBS/aligned_reads/$name.sam > /home/exx/Desktop/Jake/M32_GBS/aligned_reads/$name.bam
rm /home/exx/Desktop/Jake/M32_GBS/aligned_reads/$name.sam

done