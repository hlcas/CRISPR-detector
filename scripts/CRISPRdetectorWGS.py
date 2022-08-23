# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#	   File Name:	   CRISPRdetectorWGS.py
#	   Author:		   Lei Huang
#	   Date:		   2021.10.20
#	   E-mail:		   huanglei192@mails.ucas.ac.cn
#-------------------------------------------------
'''

import os
import sys
import time
import logging
import argparse
import textwrap

description = '''
------------------------------------------------------------------------------------------------------------------------

This script is designed to analyze whole genome sequencing (WGS) or BED format file defined targeted sequencing data.

Usage:

python CRISPRdetectorWGS.py
--sample: sample name & output directory name [required]
--e1: treatment group fq1 path [required]
--e2: treatment group fq2 path [optional]
--c1: control group fq2 path [optional]
--c2: control group fq2 path [optional]
--o: output path [default:'.']
--threads: number of threads to run sentieon minimap2 & driver module [default:1] 
--min_allele_frac: the minimum allelic fraction in treated sample [default:0.005] 
--max_fisher_pv_active: the maximum pvalue of the statistical difference between treated and untreated sample [default:0.05]
--bed: BED format file input to call variants of interested regions [optional]
--assembly: path to assembly in FASTA format : hg38.fa mm9.fa ... [required]    

------------------------------------------------------------------------------------------------------------------------
'''
parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--e1", help="treated group fq1 path",required=True)
parse.add_argument("--e2", help="treated group fq2 path",required=False)
parse.add_argument("--c1", help="control group fq1 path",required=False)
parse.add_argument("--c2", help="control group fq2 path",required=False)
parse.add_argument("--bed", help="bed format file path",required=False)
parse.add_argument("--o",help='output path',default='.',required=False)
parse.add_argument("--sample",help="sample name & output dir",required=True)
parse.add_argument("--assembly",help="genome path in fasta format",required=True)
parse.add_argument("--threads",  help="number of threads[15]",default=15,type=int)
parse.add_argument("--min_tumor_allele_frac", help="The minimum allelic fraction in treated sample",default='0.005',type=str)
parse.add_argument("--max_fisher_pv_active",help="The maximum pvalue of the statistical difference between treated and untreated sample",default='0.05',type=str)
 
args = parse.parse_args()
time0 =time.time()

e1 = os.path.abspath(args.e1)
if not os.path.exists(e1):
	sys.exit('Please check the path of treatment group fastqs.')
if args.e2 != None:
	e2 = os.path.abspath(args.e2)
	if not os.path.exists(e2):
		sys.exit('Please check the path of treatment group fastqs.')
if args.c1 != None:
	c1 = os.path.abspath(args.c1)
	if not os.path.exists(c1):
		sys.exit('Please check the path of control group fastqs.')
if args.c2 != None:
	c2 = os.path.abspath(args.c2)
	if not os.path.exists(c2):
		sys.exit('Please check the path of control group fastqs.')
		
if args.bed != None:
	interval_bed = os.path.abspath(args.bed)
	
sample_name = args.sample
fasta = os.path.abspath(args.assembly)
os.system('samtools faidx '+fasta+' && sync')

os.chdir(args.o)
os.system('mkdir -p ' + sample_name+'/temp/ && sync')
os.chdir(sample_name)

logger = logging.getLogger()
logger.setLevel(logging.INFO)
rq = time.strftime('%Y%m%d%H%M', time.localtime(time.time()))
fh = logging.FileHandler('main.log', mode='w')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
fh.setFormatter(formatter)
logger.addHandler(fh)

threads = str(args.threads)
filter_t_alt_frac = args.min_tumor_allele_frac

logger.info('Mapping treatment group fastqs to reference genome using minimap2.')
if args.e2 != None:
	os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample_name+'\\tSM:'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+e1+' '+e2+' | sentieon util sort -o temp/'+sample_name+'.tmp.bam -t '+threads+' --sam2bam -i - && sync')
else:
	os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample_name+'\\tSM:'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+e1+' | sentieon util sort -o temp/'+sample_name+'.tmp.bam -t '+threads+' --sam2bam -i - && sync')
logger.info('Finished : mapping treatment group fastqs to reference genome using minimap2.')

# Starting running minimap2 mapping reads to amplicons (control group)
if args.c1 != None:
	logger.info('Mapping control group fastqs to reference genome using minimap2.')
	if args.c2 != None:
		os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:control_'+sample_name+'\\tSM:control_'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+c1+' '+c2+' | sentieon util sort -o temp/'+sample_name+'.control.tmp.bam -t '+threads+' --sam2bam -i - && sync')
	else:
		os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:control_'+sample_name+'\\tSM:control_'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+c1+' | sentieon util sort -o temp/'+sample_name+'.control.tmp.bam -t '+threads+' --sam2bam -i - && sync')
	logger.info('Finished : mapping control group fastqs to reference genome using minimap2.')
	# Starting calling variants using sentieon driver TNscope caller
	logger.info('Calling variants.')
	param_list = ' --min_tumor_allele_frac '+filter_t_alt_frac+' --filter_t_alt_frac '+filter_t_alt_frac+' --max_fisher_pv_active '+args.max_fisher_pv_active+' --resample_depth 100000 --assemble_mode 4 --prune_factor 0'
	if args.bed == None:
		os.system('sentieon driver -t '+threads+' -r '+fasta+' -i temp/'+sample_name+'.tmp.bam -i temp/'+sample_name+'.control.tmp.bam --algo TNscope --tumor_sample '+sample_name+' --normal_sample control_'+sample_name+param_list+' variants.vcf.gz && sync')
	else:
		os.system('sentieon driver -t '+threads+' -r '+fasta+' -i temp/'+sample_name+'.tmp.bam -i temp/'+sample_name+'.control.tmp.bam --interval '+interval_bed+' --algo TNscope --tumor_sample '+sample_name+' --normal_sample control_'+sample_name+param_list+' variants.vcf.gz && sync')
else:
	logger.info('Calling variants')			
	param_list = ' --min_tumor_allele_frac '+filter_t_alt_frac+' --filter_t_alt_frac '+filter_t_alt_frac+' --resample_depth 100000 --assemble_mode 3'
	if args.bed == None:
		os.system('sentieon driver -t '+threads+' -r '+fasta+' -i temp/'+sample_name+'.tmp.bam --algo TNscope --tumor_sample '+sample_name+param_list+' temp/tmp.vcf.gz && sync')
	else:
		os.system('sentieon driver -t '+threads+' -r '+fasta+' -i temp/'+sample_name+'.tmp.bam --interval '+interval_bed+' --algo TNscope --tumor_sample '+sample_name+param_list+' temp/tmp.vcf.gz && sync')

# Filter low quality variants
raw_vcf = pd.read_csv('temp/tmp.vcf.gz',sep='\t',comment='#',header=None)
if len(raw_vcf.columns) == 10:
	raw_vcf.columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sample_name]
else:
	raw_vcf.columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sample_name,'control_'+sample_name]

SVs = raw_vcf[raw_vcf['FORMAT'] == 'GT:AD']
raw_vcf = raw_vcf[raw_vcf['FORMAT'] != 'GT:AD']
filterVCF = raw_vcf[raw_vcf['FILTER'].isin(['PASS','triallelic_site','alt_allele_in_normal'])]
SVs.to_csv(sample_name+'.sv.vcf',sep='\t',index=None)
filterVCF.to_csv(sample_name+'.snps_indels.vcf',sep='\t',index=None)

logger.info('Finished : variants called')
time1=time.time()
logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2)))
