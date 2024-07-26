#!/bin/bash

for file in /gatk/my_data/*1_R1_001_sorted_duplicates.bam; do

name="$(basename $file _R1_001_sorted_duplicates.bam)"
echo "working on $name"
gatk --java-options "-Xmx4g" HaplotypeCaller -R /gatk/my_data/Bs_v2_soft.mask.fasta -I $file -O /gatk/my_data/"$name".g.vcf -ERC GVCF

done
