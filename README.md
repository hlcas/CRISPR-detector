# CRISPR-detector
Here we propose our CRISPR-detector to facilitate the CRISPR-edited amplicon and whole genome sequencing data analysis, with functions that existing tools are not able to provide.   

CRISPR-detector brings the following four key innovations :  
1) optimized processing time allowing for hundreds of amplicons or whole genome sequencing data;   
2) integrated structural variation calling;   
3) edited and control sample co-analysis, to remove background variants not induced by gene-editing;    
4) functional and clinical consequences annotation of editing-induced mutations.  

## System requirements
### Sentieon module
Download sentieon toolkit from
https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-202010.03.tar.gz  
You may request a license by sending email to huanglei@genomics.cn

export SENTIEON_LICENSE=PATH_TO_SENTIEON/sentieon-genomics-202010.03/localhost_eval.lic  
export PATH=PATH_TO_SENTIEON/sentieon-genomics-202010.03/bin:$PATH

### Python packages
pip install biopython  
pip install pyfaidx  
pip install -U textwrap3  
conda install blast  
conda install samtools  

### ANNOVAR
Download ANNOVAR from
https://www.openbioinformatics.org/annovar/annovar_download_form.php  
perl annotate_variation.pl -downdb -webfrom annovar avdblist humandb/ -buildver hg38  
perl annotate_variation.pl -buildver hg38  -downdb -webfrom annovar refGene humandb/  
perl annotate_variation.pl -buildver hg38  -downdb -webfrom annovar clinvar_20210501 humandb/  
export PATH=PATH_TO_ANNOVAR/annovar:$PATH  
  
Organism Homo sapiens Experiment type sequencing data support variant annotations from refGene & ClinVar, other species may only support refGene annotations

## Usage
### 1. Single amplicon & pooled amplicons sequencing data analysis
python CRISPRdetectorCORE.py  
--sample, sample name & output dir  
--e1, treated group fq1 path, required = True  
--e2, treated group fq2 path, required = False  
--c1, control group fq2 path, required = False  
--c2, control group fq2 path, required = False  
--ref_fasta, single or pooled amplicon(s) sequence(s) path in fasta format, required=True    
--threads, number of threads to run sentieon minimap2 & driver module, default=1  
--anno, annotate variants with ANNOVAR [1] or not run ANNOVAR [2], required=False  
--assembly, path to assembly in fasta format : hg38.fa mm9.fa ... required=False  
--species, species : Homo_sapiens,Mus_musculus... required=False  
--db, ANNOVAR database path, required=False  

### 2. Whole genome sequencing (WGS) data analysis
python CRISPRdetectorWGS.py   
--e1, treated group fq1 path, required = True  
--e2, treated group fq2 path, required = False  
--c1, control group fq2 path, required = False  
--c2, control group fq2 path, required = False  
--sample, sample name & output dir  
--threads, number of threads to run sentieon minimap2 & driver module, default=1   
--bed, bed format file input to call variants of interest region, required=False   
--anno, annotate variants with ANNOVAR [1] or not run ANNOVAR [2], required=False  
--assembly, path to assembly in fasta format : hg38.fa mm9.fa ... required=False  
--species, species : Homo_sapiens,Mus_musculus... required=False  
--db, ANNOVAR database path, required=False  
