# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#       File Name:        CRISPRdetectorCORE.py
#       Description:      The script is designed to analyze deep-sequencing PCR products, aiming to compute CRISPR-triggered
#                         on-target efficiency.
#       Author:           Lei Huang
#       Date:             2021.10.20
#       E-mail:           huanglei192@mails.ucas.ac.cn
#-------------------------------------------------
'''
import sys
import time
import os
import textwrap
import argparse
from pyfaidx import Fasta
import pandas as pd
from Bio.Seq import Seq
#import scipy.stats as stats
import numpy as np
import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
rq = time.strftime('%Y%m%d%H%M', time.localtime(time.time()))

description = '''
------------------------------------------------------------------------------------------------------------------------

The script, supporting both paired-end and single-end reads, is designed to analyze deep-sequencing PCR products, aiming to compute CRISPR-triggered on-target efficiency.

------------------------------------------------------------------------------------------------------------------------
'''
parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--sample",help="sample name & output dir",required=True)
parse.add_argument("--e1", help="treated group fq1 path",required=True)
parse.add_argument("--e2", help="treated group fq2 path",required=False)
parse.add_argument("--c1", help="control group fq1 path",required=False)
parse.add_argument("--c2", help="control group fq2 path",required=False)
parse.add_argument("--amplicons_file", help = "Amplicons description file. This file is a tab-delimited text file with up to 5 columns (2 required)",required=True)
parse.add_argument("--threads",  help="number of threads[1]",default=1,type=int)
parse.add_argument("--anno", help="Annotate variants:1 or Not:0",default=0,type=int)
parse.add_argument("--assembly", help="assembly path in fasta format : hg38.fa mm9.fa ...",required=False)
parse.add_argument("--db", help="Annovar database path",required=False,default='/data/toolkit/annovar/')
parse.add_argument("--pvalueFilter",  help="Filter backgroud noise variants[0]",default=0,type=int)
parse.add_argument("--cleavage_offset", help='Center of quantification window to use within respect to the 3-end of the provided sgRNA sequence[-3]',default=-3,type=int)
parse.add_argument("--window_size", help="Defines the size (in bp) of the quantification window extending from the position specified by the cleavage_offset parameter in relation to the provided guide RNA sequence[0], 0 means whole amplicon analysis",default=0,type=int)
parse.add_argument("--CRISPResso", help="Running CRISPResso2[1] or not[0]",required=False,default=0,type=int)
parse.add_argument("--CPUs", help="number of CPUs[1]",default=1,type=str)
parse.add_argument("--debug", help="debug[1] or not[0]",default='0',type=str)
parse.add_argument("--sign", help="significates 0.05 or 0.01",default=0.05,type=float)
args = parse.parse_args()
# usage :

#sentieon = '/ldfssz1/ST_BIGDATA/USER/st_bigdata/Sentieon/Sentieon'

time0 =time.time()

threads = str(args.threads)
e1 = os.path.abspath(args.e1)
if args.e2 != None:
	e2 = os.path.abspath(args.e2)
sample_name = args.sample

amplicons_file_path = os.path.abspath(args.amplicons_file)

if args.c1 != None:
	c1 = os.path.abspath(args.c1)
if args.c2 != None:
	c2 = os.path.abspath(args.c2)

os.system('mkdir -p ' + sample_name+'/temp/ && sync')
amplicons_file = pd.read_csv(amplicons_file_path,sep='\t',header=None)

os.chdir(sample_name)
fh = logging.FileHandler('main.log', mode='w')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
fh.setFormatter(formatter)
logger.addHandler(fh)

def format_file(x,y,z):
	if x not in y:
		logger.error('Site '+z+' find conflicts , gRNA not in amplicon seq. Please check and make sure that gRNA and amplicon sequence in the same strand!\n')
		sys.exit('Please input the right sequence(s) and submit agian!')

if len(amplicons_file.columns) == 2:
	amplicons_file.columns = ['amplicon','amplicon_seq']
	amplicons_file['amplicon_seq'] =  amplicons_file['amplicon_seq'].apply(lambda x:x.upper())
elif len(amplicons_file.columns) == 3:
	amplicons_file.columns = ['amplicon','amplicon_seq','sgrna_seq']
	amplicons_file['amplicon_seq'] = amplicons_file['amplicon_seq'].apply(lambda x:x.upper())
	amplicons_file['sgrna_seq'] = amplicons_file['sgrna_seq'].apply(lambda x:x.upper())
	for i in range(len(amplicons_file)):
		if amplicons_file['sgrna_seq'].values[i] not in amplicons_file['amplicon_seq'].values[i]:
			sys.exit('Site '+amplicons_file['amplicon'].values[i]+' find conflicts , sgRNA not in amplicon seq. Please check!')

	with open('temp/sgRNAs.fa','w') as f:
		for i in range(len(amplicons_file)):
			f.write('>'+amplicons_file['amplicon'].values[i]+'\n'+amplicons_file['sgrna_seq'].values[i]+'\n')
	sgrna = Fasta(os.path.abspath('temp/sgRNAs.fa'))

	
with open('temp/amplicon_seq.fa','w') as f:
	for i in range(len(amplicons_file)):
		f.write('>'+amplicons_file['amplicon'].values[i]+'\n'+amplicons_file['amplicon_seq'].values[i]+'\n')

fasta = os.path.abspath('temp/amplicon_seq.fa')
amplicon_fas = Fasta(fasta)
# Starting running minimap2 mapping reads to amplicons (treated  group)
if args.e2 != None:
	#os.system(sentieon+' minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample_name+'\\tSM:'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+e1+' '+e2+' -o temp/'+sample_name+'.bam && sync')
	#os.system(sentieon+' util sort temp/'+sample_name+'.bam -o temp/'+sample_name+'.tmp.bam -t '+threads+' --sam2bam && sync')
	if args.debug == '0':
		logger.info('Starting mapping treatment group fastqs to amplicon(s) using minimap2.')
		os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample_name+'\\tSM:'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+e1+' '+e2+' | sentieon util sort -o temp/'+sample_name+'.tmp.bam -t '+threads+' --sam2bam -i - && sync')
		#os.system(sentieon + ' minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample_name+'\\tSM:'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+e1+' '+e2+' -o temp/'+sample_name+'.bam && sync')
		#os.system(sentieon + ' util sort temp/'+sample_name+'.bam -o temp/'+sample_name+'.tmp.bam -t '+threads+' --sam2bam && sync')
		logger.info('Finished : mapping treatment group fastqs to amplicon(s) using minimap2.')
else:
	if args.debug == '0':

		logger.info('Starting mapping treatment group fastq to amplicon(s) using minimap2.')
		os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample_name+'\\tSM:'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+e1+' | sentieon util sort -o temp/'+sample_name+'.tmp.bam -t '+threads+' --sam2bam -i - && sync')
		logger.info('Finished : mapping treatment group fastq to amplicon(s) using minimap2.')
window_size = args.window_size
cleavage_offset = args.cleavage_offset

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

# Starting running minimap2 mapping reads to amplicons (control group)
if args.c1 != None:
	if args.c2 != None:
		if args.debug == '0':
			logger.info('Starting mapping control group fastqs to amplicon(s) using minimap2.')
			os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:control_'+sample_name+'\\tSM:control_'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+c1+' '+c2+' | sentieon util sort -o temp/'+sample_name+'.control.tmp.bam -t '+threads+' --sam2bam -i - && sync')
			#os.system(sentieon +' minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:control_'+sample_name+'\\tSM:control_'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+c1+' '+c2+' -o temp/'+sample_name+'.control.bam && sync')
			#os.system(sentieon +' util sort temp/'+sample_name+'.control.bam -o temp/'+sample_name+'.control.tmp.bam -t '+threads+' --sam2bam && sync')
			logger.info('Finished : mapping control group fastqs to amplicon(s) using minimap2.')
		#os.system(sentieon+' minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:control_'+sample_name+'\\tSM:control_'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+c1+' '+c2+' -o temp/'+sample_name+'.control.bam && sync')
		#os.system(sentieon+' util sort temp/'+sample_name+'.control.bam -o temp/'+sample_name+'.control.tmp.bam -t '+threads+' --sam2bam && sync')

	else:
		if args.debug == '0':
			logger.info('Starting mapping control group fastq to amplicon(s) using minimap2.')
			os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:control_'+sample_name+'\\tSM:control_'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+c1+' | sentieon util sort -o temp/'+sample_name+'.control.tmp.bam -t '+threads+' --sam2bam -i - && sync')
			logger.info('Finished : mapping control group fastq to amplicon(s) using minimap2.')
	# Starting calling variants using sentieon driver TNhaplotyper2 caller
	if args.debug == '0':

		logger.info('Starting calling variants.')

		param_list = ' --min_tumor_allele_frac 0.005 --filter_t_alt_frac 0.005 --resample_depth 100000 --assemble_mode 3'
		param_list = ' --min_tumor_allele_frac 0.005 --filter_t_alt_frac 0.005 --max_fisher_pv_active 0.05 --resample_depth 100000 --assemble_mode 4 --prune_factor 0'
		
		#os.system('sentieon driver -t '+threads+' -r '+fasta+' -i temp/'+sample_name+'.tmp.bam -i temp/'+sample_name+'.control.tmp.bam --algo TNhaplotyper2 --tumor_sample '+sample_name+' --normal_sample control_'+sample_name+' --pcr_indel_model aggressive temp/tmp.vcf.gz && sync')
		os.system('sentieon driver -t '+threads+' -r '+fasta+' -i temp/'+sample_name+'.tmp.bam -i temp/'+sample_name+'.control.tmp.bam --algo TNscope --tumor_sample '+sample_name+' --normal_sample control_'+sample_name+param_list+' temp/tmp.vcf.gz && sync')
		logger.info('Finished : variants called.')
		#os.system(sentieon+' driver -t '+threads+' -r '+fasta+' -i temp/'+sample_name+'.tmp.bam -i temp/'+sample_name+'.control.tmp.bam --algo TNhaplotyper2 --tumor_sample '+sample_name+' --normal_sample control_'+sample_name+' temp/tmp.vcf.gz && sync')
else:
	if args.debug == '0':
		logger.info('Starting calling variants')
		param_list = ' --min_tumor_allele_frac 0.005 --filter_t_alt_frac 0.005 --resample_depth 100000 --assemble_mode 3 --pcr_indel_model aggressive'
		#param_list = ' '
		os.system('sentieon driver -t '+threads+' -r '+fasta+' -i temp/'+sample_name+'.tmp.bam --algo TNscope --tumor_sample '+sample_name+param_list+' temp/tmp.vcf.gz && sync')
		logger.info('Finished : variants called')

os.system('zcat temp/tmp.vcf.gz | grep -v \'##\' > temp/tmp.vcf && sync')

raw_vcf = pd.read_csv('temp/tmp.vcf',sep='\t')
raw_vcf = raw_vcf[raw_vcf['FORMAT'] != 'GT:AD']

if len(raw_vcf) != 0:	
	if args.anno == 1:
		# params runnning Annovar needed 
		assembly = args.assembly
		
		assembly_dic = {'hg38':'human','hg19':'human','mm10':'mouse','mm39':'mouse','GRCz11':'zebarfish','susScr11':'pig','Mmul_10':'monkey','bosTau9':'cow','v9.1':'frog','BDGP6.32':'fly','TAIR10':'tair','IRGSP-1.0':'rice'}
		
		annovar_db = os.path.join(args.db,assembly_dic[assembly]+'db')
		assembly_path = os.path.join(annovar_db,assembly+'.fa')
		# running blastn to get the coordinates of each amplicon on assembly fasta
		os.system('blastn -db '+assembly_path+' -query '+fasta+' -out temp/blast_out.txt -outfmt \"6 qaccver saccver sstrand sstart send\" -max_target_seqs 1 && sync')
		blast = pd.read_csv('temp/blast_out.txt',sep='\t',header=None)
		blast.columns = ['amplicon','chr','strand','start','end']
		os.system('mkdir -p temp/annovar && sync')
		# counting best-hit nums for each amplicon , if larger than 1 , annotates variants with the first match
		with open('blast_count.tab','w') as f:
			f.write('amplicon\tnum of hits\n')
			for  i in blast['amplicon'].unique():
				f.write(i+'\t'+str(len(blast[blast['amplicon'] == i]))+'\n')
	
	with open('temp/withControl.vcf','w') as g:
		with open('temp/annovar.vcf','w') as ag:
		
			g.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample_name+'\tcontrol_'+sample_name+'\n')
		
			ag.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample_name+'\t'+sample_name+'\n')
			for i in range(len(raw_vcf)):
				ampID = raw_vcf['#CHROM'].values[i]
				pos = raw_vcf['POS'].values[i]
				ref = raw_vcf['REF'].values[i]
				alt = raw_vcf['ALT'].values[i]
				g1 = raw_vcf[sample_name].values[i]
				g1AD = g1.split(':')[1]
				g1ADx = int(g1AD.split(',')[1])
				g1Rx = int(g1AD.split(',')[0])
				g1DP = str(g1ADx + g1Rx)
				gt1 = g1.split(':')[0][:3]+':'
				if len(raw_vcf.columns) == 11:
					g2 = raw_vcf['control_'+sample_name].values[i]
					g2AD = g2.split(':')[1]
					g2ADx = int(g2AD.split(',')[1])
					g2Rx = int(g2AD.split(',')[0])
					g2DP = str(g2ADx + g2Rx)
					gt2 = g2.split(':')[0][:3]+':'			
					if args.pvalueFilter == 1:
						if (g1ADx < 5) or (g2ADx < 5) or (g1Rx < 5) or (g2Rx < 5):
							pvalue = stats.fisher_exact([[g1ADx,g2ADx],[g1Rx,g2Rx]])[1]
						else:
							pvalue = stats.chi2_contingency([[g1ADx,g2ADx],[g1Rx,g2Rx]])[1]
						if pvalue < args.sign:
							g.write(ampID+'\t'+str(pos)+'\t.\t'+ref+'\t'+alt+'\t.\t.\t.\tGT:AD:DP\t'+gt1+g1AD+':'+g1DP+'\t'+gt2+g2AD+':'+g2DP+'\n')
							ag.write(ampID+'\t'+str(pos)+'\t.\t'+ref+'\t'+alt+'\t.\t.\t.\tGT:AD:DP\t'+gt1+g1AD+':'+g1DP+'\t'+gt2+g2AD+':'+g2DP+'\n')	
					else:
						g.write(ampID+'\t'+str(pos)+'\t.\t'+ref+'\t'+alt+'\t.\t.\t.\tGT:AD:DP\t'+gt1+g1AD+':'+g1DP+'\t'+gt2+g2AD+':'+g2DP+'\n')
						ag.write(ampID+'\t'+str(pos)+'\t.\t'+ref+'\t'+alt+'\t.\t.\t.\tGT:AD:DP\t'+gt1+g1AD+':'+g1DP+'\n')
				else:
					ag.write(ampID+'\t'+str(pos)+'\t.\t'+ref+'\t'+alt+'\t.\t.\t.\tGT:AD:DP\t'+gt1+g1AD+':'+g1DP+'\n')

			# pvalue < 0.05 * remained 
	#filter_vcf = pd.read_csv('temp/withControl.vcf',sep='\t')
	os.system('convert2annovar.pl -format vcf4 temp/annovar.vcf > temp/anno.avinput && sync')
	vcf_anno = pd.read_csv('temp/annovar.vcf',sep='\t')
	avinput = pd.read_csv('temp/anno.avinput',sep='\t',header=None)
	avinput.columns = ['amplicon','Start','End','Ref','Alt','het/hom','tmp','DP']
	avinput['GT:AD:DP'] = vcf_anno[sample_name]
                
	avinput['AD_Alt'] = avinput['GT:AD:DP'].apply(lambda x:int(x.split(':')[1].split(',')[1]))
                
	avinput['AF%'] = avinput['AD_Alt']*100/avinput['DP']
        
             
	#avinput = avinput.sort_values(by='AF%',ascending=False)

	def delta_len(r,a):
		if r == '-':
			return len(a)
		else:
			if a == '-':
				return 0 - len(r)
			else:
				return len(a)-len(r)
	avinput['Delta_len'] = avinput.apply(lambda row:delta_len(row['Ref'],row['Alt']),axis=1)

	if len(avinput) != 0:
		dic_window = {}
		for t in avinput['amplicon'].unique():
			avinput_rest = avinput[avinput['amplicon'] == t]
			os.system('mkdir -p '+t)	
			
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
			pos_ad_dic = {}
			pos_dp_dic = {}
			for k in range(1,len(amplicon_fas[t][:].seq)+1):
				pos_af_dic[k] = {}
				pos_ad_dic[k] = {}
				pos_dp_dic[k] = 0
				pos_af_dic[k] = {'A':0,'T':0,'C':0,'G':0,'I':0,'D':0}
				pos_ad_dic[k] = {'A':0,'T':0,'C':0,'G':0,'I':0,'D':0}

			 
			#if len(avinput_rest) != 0:
			#avinput_rest['AD_Alt'] = avinput_rest['GT:AD:DP'].apply(lambda x:int(x.split(':')[1].split(',')[1]))
			#avinput_rest['DP'] = avinput_rest['GT:AD:DP'].apply(lambda x:int(x.split(':')[2]))
			#avinput_rest['Delta_len'] = avinput_rest.apply(lambda row:delta_len(row['Ref'],row['Alt']),axis=1)
			for j in range(len(avinput_rest)):
				pos = avinput_rest['Start'].values[j] - 1
				adM = avinput_rest['AD_Alt'].values[j]
				dpM = avinput_rest['DP'].values[j]
				lenxy = avinput_rest['Delta_len'].values[j]
				freqM = avinput_rest['AF%'].values[j]
				if lenxy == 0:
					for m in avinput_rest['Alt'].values[j]:
						pos += 1
						pos_af_dic[pos][m] += freqM
						pos_ad_dic[pos][m] += adM
						pos_dp_dic[pos] = max(dpM,pos_dp_dic[pos])
				elif lenxy < 0:
					for m in avinput_rest['Ref'].values[j]:
						pos += 1
						pos_af_dic[pos]['D'] += freqM
						pos_ad_dic[pos]['D'] += adM
						pos_dp_dic[pos] = max(dpM,pos_dp_dic[pos])
				else:
					pos_af_dic[pos+1]['I'] += freqM
					pos_ad_dic[pos+1]['I'] += adM
					pos_dp_dic[pos+1] = max(dpM,pos_dp_dic[pos])

			with open(t+'/Out_insertion_deletion_substitution_locations.csv','w') as f:
				logger.info('Analysis site '+t+'.')
				if len(amplicons_file.columns) == 3:
					f.write('#sgRNA:'+str(sgRNA_start)+','+str(sgRNA_end)+'\n')
					f.write('#predicted_cleavage_position:'+str(max(sgRNA_end + cleavage_offset,1))+'\n')
					if window_size != 0:
						#f.write('#predicted_cleavage_position:'+str(max(sgRNA_end + cleavage_offset,1))+'\n')
						f.write('#window:'+str(window_start)+','+str(window_end)+'\n')
				f.write('POS,depth,substitution,deletion,insertion,indel,modified,substitution%,deletion%,insertion%,indel%,modified%\n')
				maxDP = max(pos_dp_dic.values())
				for i in range(1,len(amplicon_fas[t][:].seq)):
					if i in pos_af_dic.keys():
						vDP = pos_dp_dic[i]
						if vDP == 0:
							vDP = maxDP
						subAF = round(pos_af_dic[i]['A']+pos_af_dic[i]['G']+pos_af_dic[i]['C']+pos_af_dic[i]['T'],2)
						subAD = pos_ad_dic[i]['A']+pos_ad_dic[i]['G']+pos_ad_dic[i]['C']+pos_ad_dic[i]['T']
						delAF = round(pos_af_dic[i]['D'],2)
						delAD = pos_ad_dic[i]['D']
						insAF = round(pos_af_dic[i]['I'],2)
						insAD = pos_ad_dic[i]['I']
						indelAF = insAF+delAF
						indelAD = insAD+delAD
						modAF = indelAF+subAF
						modAD = indelAD+subAD
						f.write(str(i)+','+str(vDP)+','+str(subAD)+','+str(delAD)+','+str(insAD)+','+str(indelAD)+','+str(modAD)+','+str(subAF)+','+str(delAF)+','+str(insAF)+','+str(indelAF)+','+str(modAF)+'\n')
					else:
						f.write(str(i)+',0,0,0,0,0\n')

		
		#if len(avinput) != 0:
		avinput['AD_Alt'] = avinput['GT:AD:DP'].apply(lambda x:int(x.split(':')[1].split(',')[1]))
		avinput['AF%'] = avinput['AD_Alt']/avinput['DP']
		avinput['AF%'] = avinput['AF%'].apply(lambda x:round(x*100,2))
		avinput = avinput.sort_values(by='AF%',ascending=False)

		#if len(avinput) != 0:

		if args.anno == 1:
		
			with open('temp/lift.tmp.avinput','w') as f:
				for m in range(len(avinput)):
					blastx = blast[blast['amplicon'] == avinput['amplicon'].values[m]]
					Chr = blastx['chr'].values[0]
					sstrand = blastx['strand'].values[0]
					lift_pos = blastx['start'].values[0]
					if sstrand == 'minus':
						f.write(Chr+'\t'+str(lift_pos-avinput['End'].values[m]+1)+'\t'+str(lift_pos-avinput['Start'].values[m]+1)+'\t'+str(Seq(avinput['Ref'].values[m]).reverse_complement())+'\t'+str(Seq(avinput['Alt'].values[m]).reverse_complement())+'\t'+avinput['het/hom'].values[m]+'\t.\t'+str(avinput['DP'].values[m])+'\n')
					else:
						f.write(Chr+'\t'+str(lift_pos+avinput['Start'].values[m]-1)+'\t'+str(lift_pos+avinput['End'].values[m]-1)+'\t'+avinput['Ref'].values[m]+'\t'+avinput['Alt'].values[m]+'\t'+avinput['het/hom'].values[m]+'\t.\t'+str(avinput['DP'].values[m])+'\n')
			if 'human' in annovar_db:
				logger.info('Starting annotate variants using ANNOVAR with ClinVar and refGene database.')
				os.system('table_annovar.pl temp/lift.tmp.avinput '+annovar_db+' -buildver '+assembly+' -out temp/out -remove -protocol refGene,clinvar_20210501 -operation g,f -nastring . -csvout -polish && sync')
			else:
				logger.info('Starting annotate variants using ANNOVAR.')
				os.system('table_annovar.pl temp/lift.tmp.avinput '+annovar_db+' -buildver '+assembly+' -out temp/out -remove -protocol refGene -operation g  -nastring . -csvout -polish && sync')
	
			csvout = pd.read_csv('temp/out.'+assembly+'_multianno.csv')
			csvout = csvout.drop(['Ref','Alt','GeneDetail.refGene'],axis = 1)
			if 'control' in avinput.columns:
				avinput = avinput[['amplicon','Start','End','Ref','Alt','control','AF%']]
				avinput.columns = ['Amplicon','Amplicon_Start','Amplicon_End','Ref','Alt','Alt allele frequency (%) #Control','Alt allele frequency (%) #Treatment']
			else:
				avinput = avinput[['amplicon','Start','End','Ref','Alt','AF%']]
				avinput.columns = ['Amplicon','Amplicon_Start','Amplicon_End','Ref','Alt','Alt allele frequency (%)']
			
			outx = avinput.join(csvout)
			outx.to_csv('out_mutations.txt',sep='\t',index=None) 
		else:
			avinput = avinput[['amplicon','Start','End','Ref','Alt','AF%']]
			avinput.columns = ['Amplicon','Start','End','Ref','Alt','Alt allele frequency (%)']
			avinput.to_csv('out_mutations.txt',sep='\t',index=None)

else:
	logger.info('No variants called')
# Panel data summary 
if len(amplicon_fas.keys()) > 1:
	def modified_indel(x):
		if os.path.exists(x+'/Out_insertion_deletion_substitution_locations.csv'):
			tmp = pd.read_csv(x+'/Out_insertion_deletion_substitution_locations.csv',comment='#')
			#print(dic_window[x])
			#print(x)
			if x in dic_window.keys():
				tmp = tmp[tmp['POS'] >= dic_window[x][0]]
				tmp = tmp[tmp['POS'] <= dic_window[x][1]]
				print(str(dic_window[x][0])+'|'+str(dic_window[x][1]))
			if len(tmp) != 0:
				return str(tmp['modified%'].max())+'|'+str(tmp['indel%'].max())+'|'+str(tmp['substitution%'].max())+'|'+str(tmp['insertion%'].max())+'|'+str(tmp['deletion%'].max())+'|'+str(tmp['depth'].max())+'|'+str(tmp['modified'].max())+'|'+str(tmp['indel'].max())+'|'+str(tmp['substitution'].max())+'|'+str(tmp['insertion'].max())+'|'+str(tmp['deletion'].max())
			else:
				return '0|0|0|0|0|0|0|0|0|0|0'
		else:
			return '0|0|0|0|0|0|0|0|0|0|0'
	

	mapdf = pd.read_csv('mapping.out.tab',sep='\t')
	mapdf = mapdf.drop(len(mapdf)-1)
	mapdf = mapdf.sort_values(by='amplicon')
	mapdf.reset_index(drop=True)
	mapdf['tmpValues'] = mapdf['amplicon'].apply(modified_indel)
	mapdf['depth'] = mapdf['tmpValues'].apply(lambda x:x.split('|')[5])
	mapdf['substitution'] = mapdf['tmpValues'].apply(lambda x:x.split('|')[8])
	mapdf['insertion'] = mapdf['tmpValues'].apply(lambda x:x.split('|')[9])
	mapdf['deletion'] = mapdf['tmpValues'].apply(lambda x:x.split('|')[10])
	mapdf['indel'] = mapdf['tmpValues'].apply(lambda x:x.split('|')[7])
	mapdf['modified'] = mapdf['tmpValues'].apply(lambda x:x.split('|')[6])
	
	mapdf['substitution%'] = mapdf['tmpValues'].apply(lambda x:x.split('|')[2])
	mapdf['insertion%'] = mapdf['tmpValues'].apply(lambda x:x.split('|')[3])
	mapdf['deletion%'] = mapdf['tmpValues'].apply(lambda x:x.split('|')[4])
	mapdf['indel%'] = mapdf['tmpValues'].apply(lambda x:x.split('|')[1])
	mapdf['modified%'] = mapdf['tmpValues'].apply(lambda x:x.split('|')[0])

	mapdf = mapdf.drop('tmpValues',axis=1)
	
	def depthV(x,y):
		if x == '0':
			return y
		else:
			return x
	mapdf['depth'] = mapdf.apply(lambda row:depthV(row['depth'],row['n_reads']),axis=1)
		
	#mapdf.sort_values(by = ['indel%','depth'],inplace=True,ascending=False).to_csv('sitesComparison.txt',sep='\t',index=None)
	mapdf = mapdf.sort_values(by = ['indel%','depth'],ascending=False)
	mapdf.to_csv('sitesComparison.txt',sep='\t',index=None)

	readsQuartile = mapdf[mapdf['n_reads'] < mapdf['n_reads'].quantile(.25)]
	readsQuartile = readsQuartile[readsQuartile['n_reads'] > 0]
	reads0 = mapdf[mapdf['n_reads'] == 0]

	if len(reads0) != 0 or (len(readsQuartile) != 0):
		with open('mapped_reads_stat.txt','w') as f:
			if len(reads0) != 0:
				f.write('Site(s) with no reads mapped to : \n'+str(reads0['amplicon'].unique()).replace('\n','').replace(' ',',')+'\n')
			if len(readsQuartile) != 0:
				f.write('Site(s) with read(s) mapped but lower than 25th empirical quartile : \n'+str(readsQuartile['amplicon'].unique()).replace('\n','').replace(' ',',')+'\n')

		

os.system('zip -r '+sample_name+'.zip *txt */*csv && sync')
time1=time.time()
logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2)))

