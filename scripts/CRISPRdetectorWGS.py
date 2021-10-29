# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#       File Name:      CRISPRdetectorWGS.py
#       Description:    The script is designed to analyze whole genome sequencing data, aiming to compute call variants of input sample.
#       Author:         Lei Huang
#       Date:           2021/10/20
#       E-mail:         huanglei192@mails.ucas.ac.cn
#-------------------------------------------------
'''
import time
start=time.time()
print('Start at:')
print(time.strftime("%Y-%m-%d  %H:%M",time.localtime()))
import os
import textwrap
import argparse
import pandas as pd

description = '''
------------------------------------------------------------------------------------------------------------------------

The script, supporting both paired-end and single-end reads, is designed to analyze whole genome sequencing data, aiming to compute CRISPR-triggered on-target efficiency.

------------------------------------------------------------------------------------------------------------------------
'''
parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--sample",help="sample name & output dir",default="CRISPRdetector_O")
parse.add_argument("--e1", help="experimental group fastq 1 path",required=True)
parse.add_argument("--e2", help="experimental group fastq 2 path",required=False)
parse.add_argument("--c1", help="control group fastq 1 path",required=False)
parse.add_argument("--c2", help="control group fastq 2 path",required=False)
parse.add_argument("--assembly",help="genome path in fasta format",required=True)
parse.add_argument("--bed",help="bed file path")
parse.add_argument("--anno",help="Annotate variants:1 or Not:0",default=0,type=int)
parse.add_argument("--species",help="species : Homo_sapiens,Mus_musculus...",required=False)
parse.add_argument("--db",help="Annovar database path",required=False)
parse.add_argument("--threads",help="number of threads [1]",default=1)

args = parse.parse_args()

threads = str(args.threads)
sample_name = args.sample
fasta = args.assembly
os.system('mkdir -p ' + sample_name+'/tmp')

if args.bed != None:
	interval_bed = args.bed

tumor_bam = os.path.join(sample_name,'tmp/tumor.bam')
normal_bam = os.path.join(sample_name,'tmp/normal.bam')

# mapping reads to assembly fasta
if args.e2 != None:
	os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample_name+'\\tSM:'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+args.e1+' '+args.e2+' | sentieon util sort -o '+tumor_bam+' -t '+threads+' --sam2bam -i - && sync')
else:
	os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample_name+'\\tSM:'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+args.e1+' | sentieon util sort -o '+tumor_bam+' -t '+threads+' --sam2bam -i - && sync')
if args.c1 != None:
	if args.c2 != None:
		os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample_name+'\\tSM:'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+args.c1+' '+args.c2+' | sentieon util sort -o '+normal_bam+' -t '+threads+' --sam2bam -i - && sync')
	else:
		os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample_name+'\\tSM:'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+args.c1+' | sentieon util sort -o '+normal_bam+' -t '+threads+' --sam2bam -i - && sync')

print('Get bamfile at:')
print(time.strftime("%Y-%m-%d  %H:%M",time.localtime()))

# Running SENTIEON with or without bed
if args.c1 == None:
	if args.bed == None:
		os.system('sentieon driver -t '+threads+' -r '+fasta+' -i '+tumor_bam+' --algo TNhaplotyper2 --tumor_sample '+sample_name+' --pcr_indel_model aggressive '+os.path.join(sample_name,'tmp/tmp.vcf.gz && sync'))
	else:
		os.system('sentieon driver -t '+threads+' -r '+fasta+' -i '+tumor_bam+' --interval '+interval_bed+' --algo TNhaplotyper2 --tumor_sample '+sample_name+' '+os.path.join(sample_name,'tmp/tmp.vcf.gz')+' && sync')
else:
	if args.bed == None:
		os.system('sentieon driver -t '+threads+' -r '+fasta+' -i '+tumor_bam+' -i '+normal_bam+' --algo TNhaplotyper2 --tumor_sample '+sample_name+' --normal_sample Normal_'+sample_name+' '+os.path.join(sample_name,'tmp/tmp.vcf.gz')+' && sync')
	else:
		os.system('sentieon driver -t '+threads+' -r '+fasta+' -i '+tumor_bam+' -i '+normal_bam+' --interval '+interval_bed+' --algo TNhaplotyper2 --tumor_sample '+sample_name+' --normal_sample Normal_'+sample_name+' '+os.path.join(sample_name,'tmp/tmp.vcf.gz')+' && sync')

os.system('zcat '+os.path.join(sample_name,'tmp/tmp.vcf.gz')+'|grep -v \'##\' > '+os.path.join(sample_name,'tmp/out.vcf')+' && sync')

print('End at:')
print(time.strftime("%Y-%m-%d  %H:%M",time.localtime()))

end=time.time()
print('Running time: %s Seconds'%(round(end-start,2)))

# Annotate variants with ANNOVAR
if args.anno == 1:
	species = args.species
	assembly = fasta.split('/')[-1].replace('.fa','')
	with open('temp/out.vcf','r') as f:
		line = f.readlines()
	with open('temp/rm.out.vcf','w') as f:
		f.write(line[0])
		c = 0
		for i in range(len(vcf)):
			c += 1
			alt = vcf['ALT'].values[i]
			tmp = line[c].split('\t')
			num9 = tmp[9].split(':')
			AD = num9[1].split(',')
			if len(vcf.columns) == 10:
				c1 = 0
				for j in range(len(alt.split(','))):
					c1 += 1
					f.write(tmp[0]+'\t'+tmp[1]+'\t.\t'+tmp[3]+'\t'+alt.split(',')[j]+'\t.\t.\t.\tGT:AD:DP\t'+num9[0][:3]+':'+AD[0]+','+AD[c1]+':'+num9[3]+'\n')
			elif len(vcf.columns) == 11:
				num10 = tmp[10].split(':')
				ADx = num10[1].split(',')
				c1 = 0
				for j in range(len(alt.split(','))):
					c1 += 1
					f.write(tmp[0]+'\t'+tmp[1]+'\t.\t'+tmp[3]+'\t'+alt.split(',')[j]+'\t.\t.\t.\tGT:AD:DP\t'+num9[0][:3]+':'+AD[0]+','+AD[c1]+':'+num9[3]+'\t'+num10[1]+':'+ADx[0]+','+ADx[c1]+':'+num10[3]+'\n')

	# rm control group info to convert vcf to avinput format
	vcf = pd.read_csv('temp/rm.out.vcf',sep='\t')
	if 'Normal_'+sample_name in vcf.columns:
 		vcf.drop('Normal_'+sample_name,axis=1).to_csv('temp/anno.vcf',sep='\t',index=None)
	else:
		os.system('mv temp/rm.out.vcf temp/anno.vcf')
	
 	# convert vcf to avinput format
	os.system('convert2annovar.pl -format vcf4 '+os.path.join(sample_name,'tmp/anno.vcf')+' > '+os.path.join(sample_name,'tmp/tmp.avinput')+' && sync')
	if species == 'Homo_sapiens':
		os.system('table_annovar.pl '+os.path.join(sample_name,'tmp/tmp.avinput')+' '+args.db+' -buildver '+assembly+' -out '+os.path.join(sample_name,'out')+' -remove -protocol refGene,clinvar_20210501 -operation g,f -nastring . -csvout -polish && sync')
	else:
		os.system('table_annovar.pl '+os.path.join(sample_name,'tmp/tmp.avinput')+' '+args.db+' -buildver '+assembly+' -out '+os.path.join(sample_name,'out')+' -remove -protocol refGene -operation g  -nastring . -csvout -polish && sync')

	print('Variants annotation end at:')
	print(time.strftime("%Y-%m-%d  %H:%M",time.localtime()))
	end=time.time()
	print('Total running time: %s Seconds'%(round(end-start,2)))

