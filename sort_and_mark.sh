#!/bin/bash

for file in /gatk/my_data/*.bam; do

name="$(basename $file .bam)"
echo "working on $name"
touch /gatk/my_data/sorted_bams/"$name"_sorted.bam
gatk SortSam -I $file -O /gatk/my_data/sorted_bams/"$name"_sorted.bam -SO coordinate
gatk MarkDuplicates -I /gatk/my_data/sorted_bams/"$name"_sorted.bam -O /gatk/my_data/sorted_bams/"$name"_sorted_duplicates.bam -M /gatk/my_data/sorted_bams/duplication_metrics.txt
rm /gatk/my_data/sorted_bams/"$name"_sorted.bam
samtools index -c /gatk/my_data/sorted_bams/"$name"_sorted_duplicates.bam

done