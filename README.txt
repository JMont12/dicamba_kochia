#Written by Jake Montgomery
#started 11-1-2021
#This file outlines the commands run to conduct variant calling on a set of Bassia scoparia samples segregating for dicamba resistance (two F3 families, 4-1-1 and 4-5-10)
#the parent sequences are included as well as an F1 and F2 sample

#I started by trimmming the raw, demultiplexed reads on my laptop with trim_kochia_reads.sh this script can be found in the ~/Desktop/Jake/scripts folder
#file paths will need to be changed to run this command, but the settings for trimmomatic should be alright
~/Documents/Scripts/trim_kochia_reads.sh

#I aligned the trimmed reads to the Bs_v2_softmasked genome with BWA (genome available from https://www.biorxiv.org/content/10.1101/2023.05.26.542497v1)
bwa index ~/Desktop/Jake/M32_GBS/genome/Bs_v2_soft.mask.fasta
~/Desktop/Jake/scripts/align_M32_reads.sh
#this script aligns each sample to the reference genome, adds a read group (@RG) header line to each file with the name of the sample, and converts the sam file into a bam file.

#use docker to start up the GATK container and use a mounted volume to access data outside of the GATK container. This needs to be done as root.
#make sure that this directory includes your alignment (.bam) files and the reference genome fasta file
su
docker run -v /home/exx/Desktop/Jake/M32_GBS/aligned_reads:/gatk/my_data -it broadinstitute/gatk:4.2.0.0

#use picard tools to generate a dictionary (.dict) file for the reference
cd my_data
gatk CreateSequenceDictionary -R /gatk/my_data/Bs_v2_soft.mask.fasta -O /gatk/my_data/Bs_v2_soft.mask.dict
samtools faidx /gatk/my_data/Bs_v2_soft.mask.fasta

#call this script to sort the bams, mark duplicate reads, and make index files
mkdir sorted_bams
/gatk/my_data/sort_and_mark.sh

#remove the raw bam files and move the sorted bams with their index files to the working directory
rm /gatk/my_data/*.bam
mv gatk/my_data/sorted_bams/* gatk/my_data/

#use this shell script to go through and do the variant calling on each sample individually. I split up the samples to run them in parallel since multithreading is not yet an option in HaplotypCaller
#Each sample took ~ 20 minutes to run
nohup make_gvcf_files_0.sh > gvcf_0_out.txt 2>&1 &
nohup make_gvcf_files_1.sh > gvcf_1_out.txt 2>&1 &
nohup make_gvcf_files_2.sh > gvcf_2_out.txt 2>&1 &
nohup make_gvcf_files_3.sh > gvcf_3_out.txt 2>&1 &
nohup make_gvcf_files_4.sh > gvcf_4_out.txt 2>&1 &
nohup make_gvcf_files_5.sh > gvcf_5_out.txt 2>&1 &
nohup make_gvcf_files_6.sh > gvcf_6_out.txt 2>&1 &
nohup make_gvcf_files_7.sh > gvcf_7_out.txt 2>&1 &
nohup make_gvcf_files_8.sh > gvcf_8_out.txt 2>&1 &
nohup make_gvcf_files_9.sh > gvcf_9_out.txt 2>&1 &

#combine gvcf files in preparation for genotyping
#One argument file has all the individual vcf files with a -V flag before each one (-V M32_parent_D05_S321.g.vcf.gz -V 7710_parent_D03_S319.g.vcf.gz ....) and the other has all chromosome names with an -L flag (-L 1 -L 2 ...)
nohup gatk GenomicsDBImport --genomicsdb-workspace-path ./my_database --arguments_file ./samples.txt --arguments_file ./intervals.list -R ./Bs_v2_soft.mask.fasta > nohup_combine_gvcf.out

#run the joint genotyping
gatk GenotypeGVCFs -R Bs_v2_soft.mask.fasta -V gendb://my_database -O combined.vcf

#separate SNP and indel/mixed variant sites
gatk SelectVariants -V combined.vcf -select-type SNP -O combined_SNPs.vcf
gatk SelectVariants -V combined.vcf -select-type MIXED -O combined_mixed.vcf

#run hard filtering on the two variant types individually
#great resource for selecting flitering criteria: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
gatk VariantFiltration -V combined_SNPs.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 4.0" --filter-name "SOR4" -filter "FS > 20.0" --filter-name "FS20" -filter "MQ < 59.0" --filter-name "MQ59" -O snps_filtered_1.vcf
gatk VariantFiltration -V combined_mixed.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -O mixed_filtered.vcf

#I moved the filtered vcf files back to my laptop and used the extractGT function in vcfR to pull out just the genotypes. The code to do this is in the M32_qtl_mapping.Rmd script.
#I printed this to a file "M32_genotypes.csv"
#I used the following bash command to keep only the markers that were different and homozygous between parents
grep '77' M32_genotypes.csv > M32_homozygous_snps.csv
grep -v '77' M32_genotypes.csv| sed 's/,/\t/g' | awk '{a = $318; b=$319; if (a!=b) print $0}'| awk '{a = $318; if (a!="NA") print $0}'| awk '{a = $319; if (a!="NA") print $0}'| awk '{a = $318; if (a!="0/1") print $0}'| awk '{a = $318; if (a!="0|1") print $0}'| awk '{a = $318; if (a!="1|0") print $0}'| awk '{a = $319; if (a!="0/1") print $0}'| awk '{a = $319; if (a!="0|1") print $0}'| awk '{a = $319; if (a!="1|0") print $0}'| sed 's/\t/,/g' >> M32_homozygous_snps.csv
#I opened a head file of M32_genotypes and copied the header line with the sample names and pasted it at the top of the M32_homozygous_snps.csv file 

#I used the following command to split the samples from 4-1-1 and 4-5-10 up into their own files, each with the parents included 
cut -d ',' -f 1-212,318,319 M32_homozygous_snps.csv > M32_4-1-1_homozygous_snps.csv
cut -d ',' -f 1,214-316,318,319 M32_homozygous_snps.csv > M32_4-5-10_homozygous_snps.csv

#I ran the assign_alleles.py script on both F3 family files to translate 1's and 0's into R's and S's
~/Documents/Scripts/assign_alleles.py M32_4-1-1_homozygous_snps.csv M32_parent_S320 7710_parent_S317 4-1-1_alleles.csv
~/Documents/Scripts/assign_alleles.py M32_4-5-10_homozygous_snps.csv M32_parent_S320 7710_parent_S317 4-5-10_alleles.csv

#I trimmed off some junk at the end of the sample names from the genotyping to be congruent with the phenotyping names
sed -E 's/_[A-Z][0-9]+_R1_001//g' 4-1-1_alleles.csv
sed -E 's/_[A-Z][0-9]+_R1_001//g' 4-5-10_alleles.csv

#I copied the sample phenotype information from "Copy of Progeny test M32 dicamba resistant_07252020.xlsx" to "4-1-1_pheno.csv" and "4-5-10_pheno.csv" respectively.

#I needed to change the sample names again to be fully congruent
sed -i '' 's/7710xM32/7710-M32/g' 4-1-1_pheno.csv

#I used the following commands to generate gmap files
echo "marker,chr,pos(Mbp)" > 4-1-1_gmap.csv
cut -d ',' -f 1 4-1-1_alleles.csv| sed 's/_/,/g' | sed -E 's/([a-z]+[0-9]),([0-9]+)/\1_\2,\1,\2/g' | grep ',' >> 4-1-1_gmap.csv
echo "marker,chr,pos(Mbp)" > 4-5-10_gmap.csv
cut -d ',' -f 1 4-5-10_alleles.csv| sed 's/_/,/g' | sed -E 's/([a-z]+[0-9]),([0-9]+)/\1_\2,\1,\2/g' | grep ',' >> 4-5-10_gmap.csv

#I loaded the data into R using the 4-1-1_qtl_control_file.yaml and 4-5-10_qtl_control_file.yaml files and used the 4-1-1_qtl_mapping.Rmd script to go through a genome scan

#I reran the filtering because it looks like there is too much noise in the dataset based on the plot of genome probabilities across chromosome 4. These new filters passed ~120000 variants instead of ~160000 from the old filters
gatk VariantFiltration -V combined_SNPs.vcf -filter "QD < 5.0" --filter-name "QD5" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 2.0" --filter-name "SOR2" -filter "FS > 10.0" --filter-name "FS10" -filter "MQ < 0.0" --filter-name "MQ50" -O snps_filtered.vcf1
gatk VariantFiltration -V combined_mixed.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -O mixed_f

#I moved the new filtered snp set back to my computer and pulled out only snps passing the filters
grep '#' snps_filtered_1.vcf > snps_passing_1.vcf
grep '	PASS' snps_filtered_1.vcf >> snps_passing_1.vcf

#convert any genotypes supported by less than 5 reads to ./. This should reduce the amount of true heterozygote that get pushed to be homozygous.
~/Documents/Scripts/filter_snps_by_depth.py ./snps_passing_1.vcf 5 ./snps_passing_depth_1.vcf

#go to the M32_mapping.rmd script in Rstudio and pull out genotypes from the vcf and write them to a file called M32_genotypes_filtered.csv
#pull out only variants that are homozygous, but heterogeneous in the parents
sed -i '' 's/\"//g' M32_genotypes_filtered.csv
grep '7710-' M32_genotypes_filtered.csv > M32_homozygous_filtered_snps.csv
grep -v '7710-' M32_genotypes_filtered.csv | sed 's/,/\t/g' | awk '{a = $318; b=$319; if (a!=b) print $0}'| awk '{a = $318; if (a!="NA") print $0}'| awk '{a = $319; if (a!="NA") print $0}'| awk '{a = $318; if (a!="0/1") print $0}'| awk '{a = $318; if (a!="0|1") print $0}'| awk '{a = $318; if (a!="1|0") print $0}'| awk '{a = $319; if (a!="0/1") print $0}'| awk '{a = $319; if (a!="0|1") print $0}'| awk '{a = $319; if (a!="1|0") print $0}'| sed 's/\t/,/g' >> M32_homozygous_filtered_snps.csv

#I used the following command to split the samples from 4-1-1 and 4-5-10 up into their own files, each with the parents included
cut -d ',' -f 1-213,317-319 M32_homozygous_filtered_snps.csv > M32_4-1-1_homozygous_filtered_snps.csv
cut -d ',' -f 1,213-319 M32_homozygous_filtered_snps.csv > M32_4-5-10_homozygous_filtered_snps.csv

#I ran the assign_alleles.py script on both F3 family files to translate 1's and 0's into R's and S's
~/Documents/Scripts/assign_alleles.py M32_4-1-1_homozygous_filtered_snps.csv M32_parent_combined 7710_parent_combined 4-1-1_filtered_alleles.csv
~/Documents/Scripts/assign_alleles.py M32_4-5-10_homozygous_filtered_snps.csv M32_parent_combined 7710_parent_combined 4-5-10_filtered_alleles.csv

#I trimmed off some junk at the end of the sample names from the genotyping to be congruent with the phenotyping names
sed -i '' -E 's/_[A-Z][0-9]+_R1_001//g' 4-1-1_filtered_alleles.csv
sed -i '' -E 's/_[A-Z][0-9]+_R1_001//g' 4-5-10_filtered_alleles.csv

#you will need to change the geno file name in the control.yaml file to be 4-1-1_filtered_alleles.csv at this point
#hop back into the M32_qtl_mapping.Rmd script and try it again!

#running through this script should result in an interval scan plot similar to the one in the article. Congrats, you just mapped dicamba resistance!
