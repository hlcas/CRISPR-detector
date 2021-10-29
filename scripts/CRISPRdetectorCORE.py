# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#       File Name:        CRISPRdetectorPooled.py
#       Description:      The script is designed to analyze deep-sequencing PCR products, aiming to compute CRISPR-triggered
#                                       on-target efficiency.
#       Author:           Lei Huang
#       Date:             2021.10.20
#       E-mail:           huanglei192@mails.ucas.ac.cn
#-------------------------------------------------
'''
import time
start=time.time()
print('Start at:')
print(time.strftime("%Y-%m-%d  %H:%M",time.localtime()))
start = time.time()
import os
import textwrap
import argparse
from pyfaidx import Fasta
import pandas as pd
from Bio.Seq import Seq

description = '''
------------------------------------------------------------------------------------------------------------------------

The script, supporting both paired-end and single-end reads, is designed to analyze deep-sequencing PCR products, aiming to compute CRISPR-triggered on/off target efficiency.

------------------------------------------------------------------------------------------------------------------------
'''

parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--sample",help="sample name & output dir",required=True)
parse.add_argument("--e1", help="treated group fq1 path",required=True)
parse.add_argument("--e2", help="treated group fq2 path",required=False)
parse.add_argument("--c1", help="control group fq1 path",required=False)
parse.add_argument("--c2", help="control group fq2 path",required=False)

parse.add_argument("--ref_fasta", help="ref amplicon or hg38 ... path in fasta format",required=True)
parse.add_argument("--threads",  help="number of threads [1]",default=1,type=int)
parse.add_argument("--anno", help="Annotate variants:1 or Not:0",default=0,type=int)
parse.add_argument("--assembly", help="assembly path in fasta format : hg38.fa mm9.fa ...",required=False)
parse.add_argument("--species", help="species : Homo_sapiens,Mus_musculus...",required=False)
parse.add_argument("--db", help="Annovar database path",required=False)

args = parse.parse_args()
threads = str(args.threads)
fq1 = os.path.abspath(args.e1)

if args.e2 != None:
	fq2 = os.path.abspath(args.e2)
sample_name = args.sample
fasta = os.path.abspath(args.ref_fasta)

if args.c1 != None:
	c1 = os.path.abspath(args.c1)
if args.c2 != None:
	c2 = os.path.abspath(args.c2)

amplicon_fas = Fasta(fasta)
os.system('mkdir -p ' + sample_name+'/temp/ && sync')
os.chdir(sample_name)

# mapping reads to amplicon(s) fasta
if args.e2 != None:
	os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample_name+'\\tSM:'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+fq1+' '+fq2+' | sentieon util sort -o temp/'+sample_name+'.tmp.bam -t '+threads+' --sam2bam -i - && sync')
else:
	os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample_name+'\\tSM:'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+fq1+' | sentieon util sort -o temp/'+sample_name+'.tmp.bam -t '+threads+' --sam2bam -i - && sync')

# mapped reads nums for each amplicon
os.system('samtools idxstats temp/'+sample_name+'.tmp.bam > temp/mapping.tmp.tab && sync')

map_tab = pd.read_csv('temp/mapping.tmp.tab',sep='\t',header=None)
map_tab.columns = ['amplicons','len','mapped','unmapped']

# mapping info
with open('mapping.out.tab','w') as f:
	for i in range(len(map_tab)-1):
		f.write(map_tab['amplicons'].values[i]+'\t'+str(map_tab['mapped'].values[i])+'\n')
	#f.write('unmapped\t'+str(map_tab['unmapped'].sum())+'\n')

# mapping control group reads to amplicon(s) fasta
if args.c1 != None:
	if args.c2 != None:
		os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:control_'+sample_name+'\\tSM:control_'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+c1+' '+c2+' | sentieon util sort -o temp/'+sample_name+'.control.tmp.bam -t '+threads+' --sam2bam -i - && sync')
	else:
		os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:control_'+sample_name+'\\tSM:control_'+sample_name+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+c1+' | sentieon util sort -o temp/'+sample_name+'.control.tmp.bam -t '+threads+' --sam2bam -i - && sync')
	# sentieon driver running
	os.system('sentieon driver -r '+fasta+' -i temp/'+sample_name+'.tmp.bam -i temp/'+sample_name+'.control.tmp.bam --algo TNhaplotyper2 --tumor_sample '+sample_name+' --normal_sample control_'+sample_name+' temp/tmp.vcf.gz && sync')
else:
	os.system('sentieon driver -r '+fasta+' -i temp/'+sample_name+'.tmp.bam --algo TNhaplotyper2 --tumor_sample '+sample_name+' --pcr_indel_model aggressive temp/tmp.vcf.gz && sync')

os.system('zcat temp/tmp.vcf.gz | grep -v \'##\' > temp/tmp.vcf && sync')

vcf = pd.read_csv('temp/tmp.vcf',sep='\t')

if len(vcf) != 0:	
	if args.anno == 1:
		species = args.species
		assembly_path = args.assembly
		assembly = assembly_path.split('/')[-1].replace('.fa','')	
		# blast to compute position of amplicon sequences on assembly fasta for annotations	
		os.system('blastn -db '+assembly_path+' -query '+fasta+' -out temp/blast_out.txt -outfmt \"6 qaccver saccver sstrand sstart send\" -max_target_seqs 1 && sync')
		blast = pd.read_csv('temp/blast_out.txt',sep='\t',header=None)
		blast.columns = ['amplicon','chr','strand','start','end']
		os.system('mkdir -p temp/annovar && sync')
		# amplicon may blast to different sites on assembly fasta with best hits
		with open('blast_count.tab','w') as f:
			f.write('amplicon\tnum of hits\n')
			for  i in blast['amplicon'].unique():
				f.write(i+'\t'+str(len(blast[blast['amplicon'] == i]))+'\n')

	# vcf files alt may be split by ','
	with open('temp/tmp.vcf','r') as f:
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
	if 'control_'+sample_name in vcf.columns:
		vcf.drop('control_'+sample_name,axis=1).to_csv('temp/out.vcf',sep='\t',index=None)
	else:
		os.system('mv temp/rm.out.vcf temp/out.vcf')

	# convert vcf to avinput format
	os.system('convert2annovar.pl -format vcf4 temp/out.vcf > temp/tmp.avinput && sync')	
	avinput = pd.read_csv('temp/tmp.avinput',sep='\t',header=None)
	avinput.columns = ['amplicon','Start','End','Ref','Alt','het/hom','tmp','DP']
	vcf = pd.read_csv('temp/out.vcf',sep='\t')
	avinput['GT:AD:DP'] = vcf[sample_name]
	avinput['AD_Alt'] = avinput['GT:AD:DP'].apply(lambda x:int(x.split(':')[1].split(',')[1]))

	# indel length
	def delta_len(r,a):
		if r == '-':
			return len(a)
		else:
			if a == '-':
				return 0 - len(r)
			else:
				return len(a)-len(r)
	
	os.system('mkdir -p temp/pos/ && sync')

	for t in avinput['amplicon'].unique():
		vcfx = avinput[avinput['amplicon'] == t]
		if len(vcfx) != 0:
			os.system('mkdir '+t+' && sync') 
			vcfx['Delta_len'] = vcfx.apply(lambda row:delta_len(row['Ref'],row['Alt']),axis=1)	
			sub = vcfx[vcfx['Delta_len'] == 0]
			delx = vcfx[vcfx['Delta_len'] < 0]
			ins = vcfx[vcfx['Delta_len'] > 0]
			indel = vcfx[vcfx['Delta_len'] != 0]
			crn = vcfx['DP'].max()
			# Count indel length
			with open(t+'/Insertion_deletion_size_hist.csv','w') as f:
				f.write('Size,No.,Size%\n')
				for i in indel['Delta_len'].unique():
					nums = indel[indel['Delta_len'] == i]['AD_Alt'].sum()
					f.write(str(i)+','+str(nums)+','+str(round(nums*100/crn,2))+'\n')

			# Count AD values for each POS
			with open('temp/pos/'+t+'_POS_count.csv','w') as f:
				f.write('Type,POS,AD_Alt,DP\n')
				if len(sub) != 0:
					for i in range(len(sub)):
						if len(sub['Ref'].values[i]) == 1:
							f.write('S,'+str(sub['Start'].values[i])+','+str(sub['AD_Alt'].values[i])+','+str(sub['DP'].values[i])+'\n')
						else:
							c = -1
							for j in range(len(sub['Ref'].values[i])):
								c += 1
								f.write('S,'+str(sub['Start'].values[i]+c)+','+str(sub['AD_Alt'].values[i])+','+str(sub['DP'].values[i])+'\n')
				if len(delx) != 0:
					for i in range(len(delx)):
						c = 0
						for j in range(len(delx['Ref'].values[i])-1):
							c += 1
							f.write('D,'+str(delx['Start'].values[i]+c)+','+str(delx['AD_Alt'].values[i])+','+str(delx['DP'].values[i])+'\n')
				if len(ins) != 0:
					for i in range(len(ins)):
						f.write('I,'+str(ins['Start'].values[i])+','+str(ins['AD_Alt'].values[i])+','+str(ins['DP'].values[i])+'\n')

			# POS count table
			pos_c = pd.read_csv('temp/pos/'+t+'_POS_count.csv')

			# mutation ratio on each POS 
			with open(t+'/Out_insertion_deletion_substitution_locations.csv','w') as f:
				f.write('POS,substitution,deletion,insertion,indel,substitution%,deletion%,insertion%,indel%\n')
				for i in range(1,len(amplicon_fas[t])+1):
					if len(pos_c[pos_c['POS'] == i]) != 0:
						crn_pos = pos_c[pos_c['POS'] == i]['DP'].max()
						subn = pos_c[(pos_c['POS'] == i) & (pos_c['Type'] == 'S')]['AD_Alt'].sum()
						deln = pos_c[(pos_c['POS'] == i) & (pos_c['Type'] == 'D')]['AD_Alt'].sum()
						insn = pos_c[(pos_c['POS'] == i) & (pos_c['Type'] == 'I')]['AD_Alt'].sum()
						indel_n = insn+deln
						f.write(str(i)+','+str(subn)+','+str(deln)+','+str(insn)+','+str(indel_n)+','+str(round(subn*100/crn_pos,2))+','+str(round(deln*100/crn_pos,2))+','+str(round(insn*100/crn_pos,2))+','+str(round(indel_n*100/crn_pos,2))+'\n')
					else:
						f.write(str(i)+',0,0,0,0,0,0,0,0\n')

	if args.anno == 1:
		# generate avinput format file	
		with open('temp/lift.tmp.avinput','w') as f:
			for m in range(len(avinput)):
				blastx = blast[blast['amplicon'] == avinput['amplicon'].values[m]]
				Chr = blastx['chr'].values[0]
				sstrand = blastx['strand'].values[0]
				lift_pos = blastx['start'].values[0]
				if sstrand == 'minus':
					f.write(Chr+'\t'+str(lift_pos-avinput['Start'].values[m]+1)+'\t'+str(lift_pos-avinput['End'].values[m]+1)+'\t'+str(Seq(avinput['Ref'].values[m]))+'\t'+str(Seq(avinput['Alt'].values[m]))+'\t'+avinput['het/hom'].values[m]+'\t.\t'+str(avinput['DP'].values[m])+'\n')
				else:
					f.write(Chr+'\t'+str(lift_pos+avinput['Start'].values[m]-1)+'\t'+str(lift_pos+avinput['End'].values[m]-1)+'\t'+avinput['Ref'].values[m]+'\t'+avinput['Alt'].values[m]+'\t'+avinput['het/hom'].values[m]+'\t.\t'+str(avinput['DP'].values[m])+'\n')
		if species == 'Homo_sapiens':
			# running annovar
			os.system('table_annovar.pl temp/lift.tmp.avinput '+args.db+' -buildver '+assembly+' -out temp/out -remove -protocol refGene,clinvar_20210501 -operation g,f -nastring . -csvout -polish && sync')
		else:
			os.system('table_annovar.pl temp/lift.tmp.avinput '+args.db+' -buildver '+assembly+' -out temp/out -remove -protocol refGene -operation g  -nastring . -csvout -polish && sync')

		csvout = pd.read_csv('temp/out.'+assembly+'_multianno.csv')
		csvout = csvout.drop(['Ref','Alt'],axis = 1)
		avinput = avinput[['amplicon','Start','End','Ref','Alt','DP','AD_Alt']]
		avinput.columns = ['amplicon','ampliconStart','ampliconEnd','Ref','Alt','DP','AD_Alt']
		outx = avinput.join(csvout)
		outx.to_csv('out_mutations.tab',sep='\t',index=None) 
	else:
		# generate detailed info for each mutation
		avinput = avinput[['amplicon','Start','End','Ref','Alt','DP','AD_Alt']]
		avinput.to_csv('out_mutations.tab',sep='\t',index=None)
			
print('End at:')
print(time.strftime("%Y-%m-%d  %H:%M",time.localtime()))
end=time.time()
print('Running time: %s Seconds'%(round(end-start,2)))


