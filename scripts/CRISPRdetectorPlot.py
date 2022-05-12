# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#	   File Name:		CRISPRdetectorPlot.py
#	   Description:	  The script is designed to analyze deep-sequencing PCR products, aiming to compute CRISPR-triggered
#						 on-target efficiency.
#	   Author:		   Lei Huang
#	   Date:			 2021.10.20
#	   E-mail:		   huanglei192@mails.ucas.ac.cn
#-------------------------------------------------
'''

import os
import sys
import argparse
import textwrap
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

description = '''
------------------------------------------------------------------------------------------------------------------------

The script, supporting both paired-end and single-end reads, is designed to analyze deep-sequencing PCR products, aiming to compute CRISPR-triggered on-target efficiency.

------------------------------------------------------------------------------------------------------------------------
'''
parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--sample",help="sample name & output dir",required=True)
parse.add_argument("--o",help='output path',default='.',required=False)
parse.add_argument("--dpi",help='dpi',type=int,default=1800,required=False)
args = parse.parse_args()
sample_name = args.sample
os.chdir(args.o)
os.chdir(sample_name)

def sns_context(fontsize):
    conText={
    'axes.linewidth': 0.75,
    'grid.linewidth': 0.75,
    'lines.linewidth': 1.0,
    'lines.markersize': 3.0,
    'patch.linewidth': 1.0,
    'xtick.major.width': 0.75,
    'ytick.major.width': 0.75,
    'xtick.minor.width': 0.75,
    'ytick.minor.width': 0.75,
    'xtick.major.size': 2,
    'ytick.major.size': 2,
    'xtick.minor.size': 1.2,
    'ytick.minor.size': 1.2,
    'font.size': 7.5,
    'axes.labelsize': 8,
    'axes.titlesize': fontsize,
    'xtick.labelsize': fontsize,
    'ytick.labelsize': fontsize,
    'legend.fontsize': fontsize,
    'legend.title_fontsize': fontsize
    }
    return conText

def posterProcess(g,w,h,xlab,ylab):    
    inch_cm=2.54
    realFigHeight=w/inch_cm
    realFigWidth=h/inch_cm
    g.figure.set_size_inches(realFigHeight,realFigWidth)
    g.set_xlabel(xlabel=xlab)
    g.set_ylabel(ylabel=ylab)

sns_axes_style={
'xtick.bottom': True,
'ytick.left': True,
'axes.facecolor': '#EAEAF2',
}

##predicted_cleavage_position:113
#window:104,123

for i in os.listdir('.'):
	if os.path.exists(i+'/out_mutations_locations.csv'):
		df = pd.read_csv(i+'/out_mutations_locations.csv',comment='#',index_col=0)
		if len(df) != 0:
			if 'substitution%' in df.columns:
				df = df[['deletion%','insertion%','substitution%']]
				df.columns = ['deletion','insertion','substitution']
				maxv = max(df['substitution'].max(),df['insertion'].max(),df['deletion'].max(),1)
			else:
				df = df[['deletion%','insertion%']]
				df.columns = ['deletion','insertion']
				maxv = max(df['insertion'].max(),df['deletion'].max(),1)

			sns.set_theme(context=sns_context(7.5),style="darkgrid",font="sans-serif",palette=sns.color_palette("tab10"),rc=sns_axes_style)
			plt.ylim([0,maxv*1.1])
			plt.title('Frequency of mutations across the amplicon',fontdict={'fontsize':10})
			palette = sns.color_palette("husl", len(df.columns))
			g =sns.lineplot(data=df,palette=palette,dashes=False)

			xlab='Position'
			ylab='Frequency (%)'
			with open(i+'/out_mutations_locations.csv','r') as f:
				for j in f.readlines()[:3]:
					if j[0] == '#':
						if j[1] == 'p':
							pcp = int(j.split(': ')[1].strip())
							plt.axvline(x=pcp,ls='-',linewidth=0.8, color='purple',label='predicted cleavage position',alpha=0.5)
						elif j[1] == 'q':
							wleft = int(j.split(': ')[1].split('~')[0])
							wright = int(j.split('~')[1].strip())
							plt.axvline(x=wleft,ls='--',linewidth=0.8, color='black',label='quantification window',alpha=0.5)
							plt.axvline(x=wright,ls='--',linewidth=0.8, color='black',alpha=0.5)

			if 'pcp' in dir():
				plt.legend()

			posterProcess(g,12,8,xlab,ylab)
			g.figure.savefig(i+"/out_mutations_locations.png",dpi=args.dpi)
			plt.close()

	if os.path.exists(i+'/out_distribution_of_indel_size.csv'):
		dfis = pd.read_csv(i+'/out_distribution_of_indel_size.csv')
		if len(dfis) != 0:
			maxv = dfis['Normalized ratio%'].max()
			plt.ylim([0,maxv*1.1])
			plt.title('Indel size distribution',fontdict={'fontsize':10})
			g = sns.barplot(x=dfis['Size'],y=dfis['Normalized ratio%'],palette='autumn')
			xlab='Indel size (bp)'
			ylab='Normalized ratio (%)'
			posterProcess(g,12,8,xlab,ylab)
			g.figure.savefig(i+"/out_distribution_of_indel_size.png",dpi=args.dpi)
			plt.close()
