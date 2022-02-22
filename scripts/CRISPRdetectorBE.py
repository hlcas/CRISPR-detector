# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#	   File Name:		CRISPRdetectorBE.py
#	   Description:	  The script is designed to analyze deep-sequencing PCR products, aiming to compute CRISPR-triggered
#						 on-target efficiency.
#	   Author:		   Lei Huang
#	   Date:			 2021.10.20
#	   E-mail:		   huanglei192@mails.ucas.ac.cn
#-------------------------------------------------
'''

import os
import sys
import time
import logging
import argparse
import textwrap
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from pyfaidx import Fasta

description = '''
------------------------------------------------------------------------------------------------------------------------

The script, supporting both paired-end and single-end reads, is designed to analyze deep-sequencing PCR products, aiming to compute CRISPR-triggered on-target efficiency.

------------------------------------------------------------------------------------------------------------------------
'''
parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--e1", help="treated group fq1 path",required=True)
parse.add_argument("--e2", help="treated group fq2 path",required=False)
parse.add_argument("--c1", help="control group fq1 path",required=False)
parse.add_argument("--c2", help="control group fq2 path",required=False)
parse.add_argument("--sample",help="sample name & output dir",required=True)
parse.add_argument("--amplicons_file", help="Amplicons description file. This file is a tab-delimited text file with up to 5 columns (2 required)",required=True)
parse.add_argument("--o",help='output path',default='.',required=False)
parse.add_argument("--threads",  help="number of threads[15]",default=15,type=int)
parse.add_argument("--anno", help="Annotate variants:1 or Not:0",default=0,type=int)
parse.add_argument("--db", help="Annovar database path",required=False,default='/data/toolkit/annovar/')
parse.add_argument("--assembly", help="assembly version",required=False)
parse.add_argument("--ClinVar", help="Organism Homo sapiens experiment type sequencing data support variant annotations from ClinVar[1]", default=0, type = int)
parse.add_argument("--cleavage_offset", help='Center of quantification window to use within respect to the 3-end of the provided sgRNA sequence[-3]',default=-3,type=int)
parse.add_argument("--window_size", help="Defines the size (in bp) of the quantification window extending from the position specified by the cleavage_offset parameter in relation to the provided guide RNA sequence[0], 0 means whole amplicon analysis",default=0,type=int)
parse.add_argument("--min_tumor_allele_frac", help="The minimum allelic fraction in treated sample",default='0.005',type=str)
parse.add_argument("--min_num_of_reads",help="The minimum number of reads (per locus site) to evaluate",default=500,type=int)
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

sample_name = args.sample
amplicons_file = pd.read_csv(args.amplicons_file,sep='\t',header=None)
if len(amplicons_file.columns == 1):
	amplicons_file = pd.read_csv(args.amplicons_file,sep='\s+',header=None)
	if len(amplicons_file.columns) == 1 or (len(amplicons_file.columns) > 3):
		sys.exit('Please check the format of input amplicons description file.')
elif len(amplicons_file.columns) > 3:
	sys.exit('Please check the format of input amplicons description file.')

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
window_size = args.window_size
cleavage_offset = args.cleavage_offset
min_num_of_reads = args.min_num_of_reads
filter_t_alt_frac = args.min_tumor_allele_frac

def format_file(x,y,z):
	if x not in y:
		if x not in str(Seq(y).reverse_complement()):
			logger.error('Site '+z+' find conflicts , gRNA not in amplicon seq!\n')
			sys.exit('Please input the right sequence(s) and submit agian!')
		else:
			return str(Seq(y).reverse_complement())
	else:
		return y

if len(amplicons_file.columns) == 2:
	amplicons_file.columns = ['amplicon','amplicon_seq']
	amplicons_file['amplicon_seq'] =  amplicons_file['amplicon_seq'].apply(lambda x:x.upper())
elif len(amplicons_file.columns) == 3:
	amplicons_file.columns = ['amplicon','amplicon_seq','sgrna_seq']
	amplicons_file['amplicon_seq'] = amplicons_file['amplicon_seq'].apply(lambda x:x.upper())
	amplicons_file['sgrna_seq'] = amplicons_file['sgrna_seq'].apply(lambda x:x.upper())
	amplicons_file['amplicon_seq'] = amplicons_file.apply(lambda row:format_file(row['sgrna_seq'],row['amplicon_seq'],row['amplicon']),axis=1)

	with open('temp/sgRNAs.fa','w') as f:
		for i in range(len(amplicons_file)):
			f.write('>'+amplicons_file['amplicon'].values[i]+'\n'+amplicons_file['sgrna_seq'].values[i]+'\n')
	sgrna = Fasta(os.path.abspath('temp/sgRNAs.fa'))

with open('temp/amplicon_seq.fa','w') as f:
	for i in range(len(amplicons_file)):
		f.write('>'+amplicons_file['amplicon'].values[i]+'\n'+amplicons_file['amplicon_seq'].values[i]+'\n')

fasta = os.path.abspath('temp/amplicon_seq.fa')
amplicon_fas = Fasta(fasta)

logger.info('Mapping treatment group fastqs to amplicon(s) using minimap2.')
if args.e2 != None:
	os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample_name+'\\tSM:'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+e1+' '+e2+' | sentieon util sort -o temp/'+sample_name+'.tmp.bam -t '+threads+' --sam2bam -i - && sync')
else:
	os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample_name+'\\tSM:'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+e1+' | sentieon util sort -o temp/'+sample_name+'.tmp.bam -t '+threads+' --sam2bam -i - && sync')
logger.info('Finished : mapping treatment group fastqs to amplicon(s) using minimap2.')

# Numbers of reads mapped to each amplicon
os.system('samtools idxstats temp/'+sample_name+'.tmp.bam > temp/mapping.tmp.tab && sync')

map_tab = pd.read_csv('temp/mapping.tmp.tab',sep='\t',header=None)
map_tab.columns = ['amplicons','len','mapped','unmapped']

# mapping info
with open('mapping.out.tab','w') as f:
	f.write('amplicon\tn_reads\n')
	for i in range(len(map_tab)-1):
		f.write(map_tab['amplicons'].values[i]+'\t'+str(map_tab['mapped'].values[i])+'\n')
	f.write('unmapped\t'+str(map_tab['unmapped'].sum())+'\n')

mapdf = pd.read_csv('mapping.out.tab',sep='\t')
mapdf = mapdf.drop(len(mapdf)-1)

# Starting running minimap2 mapping reads to amplicons (control group)
if args.c1 != None:
	logger.info('Mapping control group fastqs to amplicon(s) using minimap2.')
	if args.c2 != None:
		os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:control_'+sample_name+'\\tSM:control_'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+c1+' '+c2+' | sentieon util sort -o temp/'+sample_name+'.control.tmp.bam -t '+threads+' --sam2bam -i - && sync')
	else:
		os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:control_'+sample_name+'\\tSM:control_'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+c1+' | sentieon util sort -o temp/'+sample_name+'.control.tmp.bam -t '+threads+' --sam2bam -i - && sync')
	os.system('samtools idxstats temp/'+sample_name+'.control.tmp.bam > temp/control.mapping.tab && sync')
	mapdf_c = pd.read_csv('temp/control.mapping.tab',sep='\t',header=None)
	mapdf_c.columns = ['amplicon','len','n_reads|control','unmapped']
	mapdf_c = mapdf_c[['amplicon','n_reads|control']]
	mapdf = pd.merge(mapdf,mapdf_c,on='amplicon',how='left')
	logger.info('Finished : mapping control group fastqs to amplicon(s) using minimap2.')
	# Starting calling variants using sentieon driver TNscope caller
	logger.info('Calling variants.')
	param_list = ' --keep_overlap 1 --min_tumor_allele_frac '+filter_t_alt_frac+' --filter_t_alt_frac '+filter_t_alt_frac+' --max_fisher_pv_active '+args.max_fisher_pv_active+' --resample_depth 100000 --assemble_mode 4 --prune_factor 0'
	os.system('sentieon driver -t '+threads+' -r '+fasta+' -i temp/'+sample_name+'.tmp.bam -i temp/'+sample_name+'.control.tmp.bam --algo TNscope --tumor_sample '+sample_name+' --normal_sample control_'+sample_name+param_list+' temp/tmp.vcf.gz && sync')
else:
	logger.info('Calling variants')			
	param_list = ' --keep_overlap 1 --min_tumor_allele_frac '+filter_t_alt_frac+' --filter_t_alt_frac '+filter_t_alt_frac+' --resample_depth 100000 --assemble_mode 3'
	os.system('sentieon driver -t '+threads+' -r '+fasta+' -i temp/'+sample_name+'.tmp.bam --algo TNscope --tumor_sample '+sample_name+param_list+' temp/tmp.vcf.gz && sync')

logger.info('Finished : variants called')

try:
	raw_vcf = pd.read_csv('temp/tmp.vcf.gz',sep='\t',header=None,comment='#',compression='gzip')
except:
	logger.info('No variants called.')
	time1=time.time()
	logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2)))
	sys.exit('No variants called.')

if len(raw_vcf.columns) == 10:
	raw_vcf.columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sample_name]
else:
	raw_vcf.columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sample_name,'control_'+sample_name]

raw_vcf['tGT'] = raw_vcf[sample_name].apply(lambda x:x.split(':')[0]+':')
raw_vcf['tAF'] = raw_vcf[sample_name].apply(lambda x:x.split(':')[2]+':')
raw_vcf['tDP'] = raw_vcf[sample_name].apply(lambda x:int(x.split(':')[3]))

def clean_reads(x,y,z):
	if len(raw_vcf[raw_vcf['#CHROM'] == x]) == 0:
		return y
	else:
		return raw_vcf[raw_vcf['#CHROM'] == x][z].max()

# clean reads nums Mock number of reads	Treatment number of reads

mapdf['n_reads'] = mapdf.apply(lambda row:clean_reads(row['amplicon'],row['n_reads'],'tDP'),axis=1)
loci_low_reads = list(mapdf[mapdf['n_reads'] < min_num_of_reads]['amplicon'].values)
mapdf = mapdf[mapdf['n_reads'] >= min_num_of_reads]
loci_high_reads = mapdf['amplicon'].values
mapdf = mapdf.reset_index(drop=True)

def vcflencheck(vcfx):
	if len(vcfx) == 0:
		logger.info('No variants called.')
		time1=time.time()
		logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2)))
		sys.exit('No variants called.')
	
# filter out sites with clean reads num < param min_num_of_reads
raw_vcf = raw_vcf[raw_vcf['#CHROM'].isin(loci_high_reads)]
vcflencheck(raw_vcf)
raw_vcf = raw_vcf[raw_vcf['FORMAT'] != 'GT:AD']
vcflencheck(raw_vcf)
sv_vcf	= raw_vcf[raw_vcf['FORMAT'] == 'GT:AD']

def differ_len(r,a):
	if r == '-':
		return len(a)
	else:
		if a == '-':
			return 0 - len(r)
		else:
			return len(a)-len(r)


raw_vcf['differ_len'] = raw_vcf.apply(lambda row:differ_len(row['REF'],row['ALT']),axis=1)
if args.ignore_substitutions == 1:
	raw_vcf = raw_vcf[raw_vcf['differ_len'] != 0]
	vcflencheck(raw_vcf)

raw_vcf = raw_vcf.reset_index(drop=True)

with open('temp/treatment.vcf','w') as f:		
	f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample_name+'\n')
	for i in range(len(raw_vcf)):
		f.write(raw_vcf['#CHROM'].values[i]+'\t'+str(raw_vcf['POS'].values[i])+'\t.\t'+raw_vcf['REF'].values[i]+'\t'+raw_vcf['ALT'].values[i]+'\t.\t.\t.\tGT:AF:DP\t'+raw_vcf['tGT'].values[i]+raw_vcf['tAF'].values[i]+str(raw_vcf['tDP'].values[i])+'\n')

if args.c1 != None:		
	raw_vcf['cGT'] = raw_vcf['control_'+sample_name].apply(lambda x:x.split(':')[0]+':')
	raw_vcf['cAF'] = raw_vcf['control_'+sample_name].apply(lambda x:x.split(':')[2]+':')
	raw_vcf['cDP'] = raw_vcf['control_'+sample_name].apply(lambda x:x.split(':')[3])

	with open('temp/treatment.vcf','r') as f:
		lines = f.readlines()
	with open('temp/paired.vcf','w') as g:
		g.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample_name+'\tcontrol_'+sample_name+'\n')
		c = 0
		for i in range(len(raw_vcf)):		
			c += 1
			g.write(lines[c].strip()+'\t'+raw_vcf['cGT'].values[i]+raw_vcf['cAF'].values[i]+raw_vcf['cDP'].values[i]+'\n')
	mapdf['n_reads|control'] = mapdf.apply(lambda row:clean_reads(row['amplicon'],row['n_reads|control'],'cDP'),axis=1)
	mapdf.columns = ['n_reads|treatment' if i == 'n_reads' else i for i in mapdf.columns]

os.system('convert2annovar.pl -format vcf4 temp/treatment.vcf > temp/anno.avinput && sync')
vcf_anno = pd.read_csv('temp/treatment.vcf',sep='\t')
avinput = pd.read_csv('temp/anno.avinput',sep='\t',header=None)
avinput.columns = ['amplicon','start','end','ref','alt','het/hom','tmp','DP']
avinput['GT:AF:DP'] = vcf_anno[sample_name]			  
avinput['AF%'] = avinput['GT:AF:DP'].apply(lambda x:float(x.split(':')[1])*100)

if os.path.exists('temp/paired.vcf'):
	avinput['control'] = raw_vcf['cAF'].apply(lambda x:round(float(x[:-1])*100,2))

avinput['differ_len'] = avinput.apply(lambda row:differ_len(row['ref'],row['alt']),axis=1)

if len(avinput) != 0:
	dic_window = {}

	def indelSIZE(vx,vy,vz):
		if vx not in vy:
			vy[vx] = vz
		else:
			vy[vx] += vz

	for t in avinput['amplicon'].unique():
		avinput_loci = avinput[avinput['amplicon'] == t]
		os.system('mkdir -p '+t)	
		avinput_loci = avinput_loci.reset_index(drop=True)	
		
		if len(amplicons_file.columns) == 3:
			sgRNA_start = amplicon_fas[t][:].seq.find(sgrna[t][:].seq) + 1
			sgRNA_end = sgRNA_start + len(sgrna[t][:].seq) - 1
			if window_size != 0:
				cleavage_offset = args.cleavage_offset
				window_start = max(sgRNA_end - window_size + cleavage_offset,1)
				window_end = min(window_start + 2*window_size,len(amplicon_fas[t][:].seq))	
			else:
				window_start = 1
				window_end = len(amplicon_fas[t][:].seq)
		else:
			window_start = 1
			window_end = len(amplicon_fas[t][:].seq)

		dic_window[t] = [window_start,window_end]

		pos_af_dic = {}
		for k in range(1,len(amplicon_fas[t][:].seq)+1):
			pos_af_dic[k] = {}
			pos_af_dic[k] = {'A':0,'T':0,'C':0,'G':0,'I':0,'D':0}

		indelsize = {}
		for j in range(len(avinput_loci)):
			pos = avinput_loci['start'].values[j] - 1
			dlen = avinput_loci['differ_len'].values[j]
			afx = avinput_loci['AF%'].values[j]

			if dlen == 0:
				for m in avinput_loci['alt'].values[j]:
					pos += 1
					pos_af_dic[pos][m] += afx
			elif dlen < 0:
				if window_size != 0:
					if len(set(range(pos + 1,pos+1-dlen)).intersection(set(range(window_start,window_end+1)))) != 0:
						indelSIZE(dlen,indelsize,afx)
				else:
					indelSIZE(dlen,indelsize,afx)

				for m in avinput_loci['ref'].values[j]:
					pos += 1
					pos_af_dic[pos]['D'] += afx
			else:
				if pos + 1 in range(window_start,window_end+1):
					indelSIZE(dlen,indelsize,afx)
				pos_af_dic[pos+1]['I'] += afx

		with open(t+'/out_mutations_locations.csv','w') as f:
			with open(t+'/out_be_mutations_locations.csv','w') as fb:

				logger.info('Analysis site '+t+'.')
				if len(amplicons_file.columns) == 3:
					f.write('#sgRNA:'+str(sgRNA_start)+','+str(sgRNA_end)+'\n')
					fb.write('#sgRNA:'+str(sgRNA_start)+','+str(sgRNA_end)+'\n')
					f.write('#predicted_cleavage_position:'+str(max(sgRNA_end + cleavage_offset,1))+'\n')
					if window_size != 0:
						f.write('#window:'+str(window_start)+','+str(window_end)+'\n')
						fb.write('#window:'+str(window_start)+','+str(window_end)+'\n')
				f.write('POS,substitution%,deletion%,insertion%,indel%,modified%\n')
				fb.write('POS,Nucleotide,A,G,C,T,-\n')

				for i in range(1,len(amplicon_fas[t][:].seq)):
					if i in pos_af_dic.keys():
						afA = round(pos_af_dic[i]['A'],2)
						afG = round(pos_af_dic[i]['G'],2)
						afC = round(pos_af_dic[i]['C'],2)
						afT = round(pos_af_dic[i]['T'],2)
						subAF = round(pos_af_dic[i]['A']+pos_af_dic[i]['G']+pos_af_dic[i]['C']+pos_af_dic[i]['T'],2)
						delAF = round(pos_af_dic[i]['D'],2)
						insAF = round(pos_af_dic[i]['I'],2)
						indelAF = insAF+delAF
						modAF = indelAF+subAF
						
						f.write(str(i)+','+str(subAF)+','+str(delAF)+','+str(insAF)+','+str(indelAF)+','+str(modAF)+'\n')
						nuclx = str(amplicon_fas[t][i-1])
						nuclxAF = round(100-subAF-delAF,2)
						if  nuclx == 'A':
							fb.write(str(i)+','+nuclx+','+str(nuclxAF)+','+str(afG)+','+str(afC)+','+str(afT)+','+str(delAF)+'\n')
						elif nuclx == 'G':
							fb.write(str(i)+','+nuclx+','+str(afA)+','+str(nuclxAF)+','+str(afC)+','+str(afT)+','+str(delAF)+'\n')
						elif nuclx == 'C':
							fb.write(str(i)+','+nuclx+','+str(afA)+','+str(afG)+','+str(nuclxAF)+','+str(afT)+','+str(delAF)+'\n')
						else:
							fb.write(str(i)+','+nuclx+','+str(afA)+','+str(afG)+','+str(afC)+','+str(nuclxAF)+','+str(delAF)+'\n')

		with open(t+'/out_distribution_of_indel_size.csv','w') as fi:
			indel_af_sum = sum(indelsize.values())
			fi.write('Size,Normalized ratio%\n')
			for ki in sorted(indelsize):
				fi.write(str(ki)+','+str(round(indelsize[ki]*100/indel_af_sum,2))+'\n')
				
	avinput['AF%'] = avinput['AF%'].apply(lambda x:round(x,2))
	avinput = avinput.sort_values(by='AF%',ascending=False)

	if args.anno == 1:

		# params runnning Annovar needed 
		assembly = args.assembly
		assembly_path = os.path.join(annovar_db,assembly+'.fa')
		# running blastn to get the coordinates of each amplicon on assembly fasta
		os.system('blastn -db '+assembly_path+' -query '+fasta+' -out temp/blast_out.txt -outfmt \"6 qaccver saccver sstrand sstart send\" -max_target_seqs 1 && sync')
		blast = pd.read_csv('temp/blast_out.txt',sep='\t',header=None)
		blast.columns = ['amplicon','chr','strand','start','end']
	
	# counting best-hit nums for each amplicon , if larger than 1 , annotates variants with the first match
		with open('blast_count.tab','w') as f:
			f.write('amplicon\tnum of hits\n')
			for  i in blast['amplicon'].unique():
				f.write(i+'\t'+str(len(blast[blast['amplicon'] == i]))+'\n')
		with open('temp/lift.tmp.avinput','w') as f:
			for m in range(len(avinput)):
				blastx = blast[blast['amplicon'] == avinput['amplicon'].values[m]]
				Chr = blastx['chr'].values[0]
				sstrand = blastx['strand'].values[0]
				lift_pos = blastx['start'].values[0]
				if sstrand == 'minus':
					f.write(Chr+'\t'+str(lift_pos-avinput['end'].values[m]+1)+'\t'+str(lift_pos-avinput['start'].values[m]+1)+'\t'+str(Seq(avinput['ref'].values[m]).reverse_complement())+'\t'+str(Seq(avinput['alt'].values[m]).reverse_complement())+'\t'+avinput['het/hom'].values[m]+'\t.\t'+str(avinput['DP'].values[m])+'\n')
				else:
					f.write(Chr+'\t'+str(lift_pos+avinput['start'].values[m]-1)+'\t'+str(lift_pos+avinput['end'].values[m]-1)+'\t'+avinput['ref'].values[m]+'\t'+avinput['alt'].values[m]+'\t'+avinput['het/hom'].values[m]+'\t.\t'+str(avinput['DP'].values[m])+'\n')

		if args.ClinVar ==1:
			logger.info('Starting annotate variants using ANNOVAR with ClinVar and refGene database.')
			os.system('table_annovar.pl temp/lift.tmp.avinput '+annovar_db+' -buildver '+assembly+' -out temp/out -remove -protocol refGene,clinvar_20210501 -operation g,f -nastring . -csvout -polish && sync')
		else:
			logger.info('Starting annotate variants using ANNOVAR.')
			os.system('table_annovar.pl temp/lift.tmp.avinput '+annovar_db+' -buildver '+assembly+' -out temp/out -remove -protocol refGene -operation g  -nastring . -csvout -polish && sync')
	
		csvout = pd.read_csv('temp/out.'+assembly+'_multianno.csv')
		csvout = csvout.drop(['Ref','Alt','GeneDetail.refGene'],axis = 1)
		#avinput = avinput.join(csvout)
	if 'control' in avinput.columns:                
		avinput = avinput[['amplicon','start','end','ref','alt','control','AF%']]
		avinput['control'] = avinput['control'].apply(lambda x:round(x,2))
		avinput.columns = ['Amplicon','Amplicon_Start','Amplicon_End','Ref','Alt','Alt allele frequency (%) | Control','Alt allele frequency (%) | Treatment']
	else:
		avinput = avinput[['amplicon','start','end','ref','alt','AF%']]
		avinput.columns = ['Amplicon','Amplicon_Start','Amplicon_End','Ref','Alt','Alt allele frequency (%)']

	if 'csvout' in dir():
		avinput = avinput.join(csvout)
	else:
		avinput.columns = [i.replace('Amplicon_','') for i in avinput.columns]

	# Filter out mutations (not overlap quantification window)
	def overlap(x,y,z):
		if len(set(range(x,y+1)).intersection(set(range(dic_window[z][0],dic_window[z][1]+1)))) != 0:
			return 'T'
		else:
			return 'F'

	for amp in avinput['Amplicon'].unique():		
		avinput_amp = avinput[avinput['Amplicon'] == amp]
		if window_size != 0:
			avinput_amp['Overlap'] = avinput_amp.apply(lambda row:overlap(row['Start'],row['End'],amp),axis=1)
			avinput_amp = avinput_amp[avinput_amp['Overlap'] == 'T']
			avinput_amp.drop('Overlap',axis=1,inplace=True)
		avinput_amp.to_csv(amp+'/out_mutations.txt',sep='\t',index=None)

# SV
if len(sv_vcf) != 0:
	sv_vcf = sv_vcf[sv_vcf['#CHROM'].isin(loci_high_reads)]
	if len(sv_vcf) != 0:
		# sg1-ON	68	ins_0	T	<INS>	47	PASS	SVTYPE=INS;SOMATIC	GT:AD	./.:4873960,672	./.:.
		sv_vcf = sv_vcf[['#CHROM','POS','ID','REF','ALT',sample_name]]
		sv_vcf['AD']  = sv_vcf[sample_name].apply(lambda x:x.split(':')[1])
		sv_vcf['ALT_AD'] = sv_vcf['AD'].apply(lambda x:int(x.split(',')[1]))
		sv_vcf['REF_AD'] = sv_vcf['AD'].apply(lambda x:int(x.split(',')[0]))
		sv_vcf['DP'] = sv_vcf['ALT_AD'] + sv_vcf['REF_AD']
		sv_vcf['AF%'] = sv_vcf['ALT_AD']/sv_vcf['REF_AD']
		sv_vcf = sv_vcf[sv_vcf['AF%'] >= filter_t_alt_frac]
		if len(sv_vcf) != 0:
			sv_vcf['AF%'] = sv_vcf['AF%'].apply(lambda x:round(x*100,2))
			#CHROM,POS,ID,REF,ALT,AF%
			sv_vcf = sv_vcf.drop([sample_name,'AD','ALT_AD','REF_AD','DP'],axis=1)
			sv_vcf.columns = ['Amplicon','Pos','SV ID','Ref','Alt','Alt allele frequency (%)']
			for i in sv_vcf['Amplicon'].unique():
				if not os.path.exists(i):
					os.system('mkdir '+i)
				sv_vcf[sv_vcf['Amplicon'] == i].to_csv('out_sv.txt',sep='\t',index=None)

# Panel data summary 
if len(amplicon_fas.keys()) > 0:
	def substitutions_removal(xloci,ydf):
		if args.ignore_substitutions == 1:
			indeldf = ydf[['POS','deletion%','insertion%','indel%']]		
			os.system('grep \'#\' '+xloci+'/out_mutations_locations.csv > '+xloci+'/tmp1.txt')
			indeldf.to_csv(xloci+'/tmp2.txt',index=None)
			os.system('cat '+xloci+'/tmp* > '+xloci+'/out_mutations_locations.csv')
			os.system('rm -rf '+xloci+'/tmp*')
		
	def modified_indel(x):
		if os.path.exists(x+'/out_mutations_locations.csv'):
			tmp = pd.read_csv(x+'/out_mutations_locations.csv',comment='#')
			if x in dic_window.keys():
				tmp = tmp[tmp['POS'] >= dic_window[x][0]]
				tmp = tmp[tmp['POS'] <= dic_window[x][1]]
				print(str(dic_window[x][0])+'|'+str(dic_window[x][1]))
			if len(tmp) != 0:
				substitutions_removal(x,tmp)
				return str(tmp['modified%'].max())+'|'+str(tmp['indel%'].max())+'|'+str(tmp['substitution%'].max())+'|'+str(tmp['insertion%'].max())+'|'+str(tmp['deletion%'].max())
			else:
				substitutions_removal(x,tmp)
				return '0|0|0|0|0'
		else:
			return '0|0|0|0|0'
		
	mapdf['tmpValues'] = mapdf['amplicon'].apply(modified_indel)
	mapdf['substitution%'] = mapdf['tmpValues'].apply(lambda x:x.split('|')[2])
	mapdf['insertion%'] = mapdf['tmpValues'].apply(lambda x:x.split('|')[3])
	mapdf['deletion%'] = mapdf['tmpValues'].apply(lambda x:x.split('|')[4])
	mapdf['indel%'] = mapdf['tmpValues'].apply(lambda x:x.split('|')[1])
	mapdf['modified%'] = mapdf['tmpValues'].apply(lambda x:x.split('|')[0])
	mapdf = mapdf.drop('tmpValues',axis=1)

	if args.ignore_substitutions == 1:
		mapdf = mapdf.drop(['substitution%','modified%'],axis=1)

	mapdf.to_csv('sitesComparison.txt',sep='\t',index=None)

	if len(loci_low_reads) != 0:
		with open('mapped_reads_stat.txt','w') as f:
			f.write('Site(s) with less than '+str(min_num_of_reads)+' reads mapped to : '+str(loci_low_reads).replace('[','').replace(']','').replace('\'','')+'.\n')
 
os.system('zip -r '+sample_name+'.zip *txt */*csv */*txt && sync')
time1=time.time()
logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2)))
