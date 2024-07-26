#!/usr/bin/python

#this script will take in a .vcf file and convert any genotype calls that are supported by less reads than some specified threshold to ./.
#this script was written to work with VCFv4.2 format output by VariantFiltration tool within GATK. I will try to keep up with changes in format convention
#please post questions/issues either to my github page or email me directly @ jake.montgomery@colostate.edu
#usage /path/to/filter_snps_by_depth.py /data/in.vcf min_depth /path/to/out.vcf

from sys import argv, version_info
from os.path import realpath, splitext
import os

in_vcf, min_depth, out_vcf = realpath(argv[1]), argv[2], realpath(argv[3])

in_fh=open(in_vcf, 'r')
out_fh=open(out_vcf, 'w+')
parts=[]
subparts=[]
col_count=0
sub_count=0
allele_depths=[]
line_count=0

#loop through the input vcf file
for line in in_fh:
	line=line.strip('\n')
	col_count=0
	if line_count==0:
		out_fh.write(str(line))
	elif line.startswith('#'):
		out_fh.write('\n'+str(line))
	else:
		out_fh.write('\n')
		parts=line.split('\t')
		for sample in parts:
			sub_count=0
			if col_count==0:
				out_fh.write(str(sample))
			elif col_count<9:
				out_fh.write('\t'+str(sample))
			else:
				subparts=sample.split(':')
				allele_depths=subparts[1].split(',')
				if int(allele_depths[0])+int(allele_depths[1])<int(min_depth):
					subparts[0]='./.'
				out_fh.write('\t')
				sub_count=0
				for section in subparts:
					if sub_count==0:
						out_fh.write(subparts[sub_count])
					else:
						out_fh.write(':'+str(subparts[sub_count]))
					sub_count+=1
			col_count+=1
	line_count+=1