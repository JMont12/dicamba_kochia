#!/bin/bash

for file in /Users/jake/Documents/PhD/Dicamba_kochia/11-1-21_GBS/Gaines_Project_002/raw_reads/*.fastq.gz; do

name="$(basename $file)"
echo "working on $name"
touch /Users/jake/Documents/PhD/Dicamba_kochia/11-1-21_GBS/Gaines_Project_002/trimmed_reads/$name
java -jar ~/apps/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 7 $file /Users/jake/Documents/PhD/Dicamba_kochia/11-1-21_GBS/Gaines_Project_002/trimmed_reads/$name ILLUMINACLIP:/Users/jake/apps/Trimmomatic-0.36/adapters/TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

done
