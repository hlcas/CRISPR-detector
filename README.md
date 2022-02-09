# CRISPR-detector

CRISPR-detector provides a web-hosted platform (https://db.cngb.org/crispr-detector/) and local deployable pipeline to fast and accurately identify and annotate editing-induced mutations from genome editing assays. 

CRISPR-detector pipeline possesses 4 key innovations :  

1) optimized processing speed allowing for hundreds of amplicons and whole genome scanning;  

2) edited and control sample co-analysis, to remove background variants exist prior to genome editing;  

3) functional and clinical annotation of editing-induced mutations;  

4) integrated structural variation (SV) calling.  


## System requirements
### Sentieon module
Download sentieon toolkit from
https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-202010.03.tar.gz  
You may request a license by sending emails to huanglei@genomics.cn

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
  
Organism Homo sapiens experiment type sequencing data support variant annotations from refGene & ClinVar, other species may only support refGene annotations

#### You may build ANNOVAR database youself for any species with coresponding genome assembly and gff3 format files
For example, to build a database for zebrafish. Download GRCz11.fa and GRCz11.gff3 from public database.  
Then running commands as following:  

conda install -c bioconda/label/cf201901 gffread  
conda install -c bioconda/label/cf201901 ucsc-gtftogenepred  
conda install -c bioconda/label/cf201901 blast  

cd PATH_TO_ANNOVAR/  
mkdir zebrafishdb && cd zebrafishdb  
mv */GRCz11.fa zebrafishdb  
mv */GRCz11.gff3 zebrafishdb  

gffread GRCz11.gff3 -T -o GRCz11.gtf  
gtfToGenePred -genePredExt GRCz11.gtf GRCz11_refGene.txt  
retrieve_seq_from_fasta.pl --format refGene --seqfile GRCz11.fa GRCz11_refGene.txt --out GRCz11_refGeneMrna.fa    
makeblastdb -in GRCz11.fa -dbtype nucl  

## Usage
### 1. Single amplicon & pooled amplicons sequencing data analysis
python scripts/CRISPRdetectorCORE.py  

--sample: sample name & output dir  

--e1: treatment group fq1 path, required = True  

--e2: treatment group fq2 path, required = False  

--c1: control group fq2 path, required = False  

--c2: control group fq2 path, required = False  

--sample: sample name & output dir name, required = True  

--amplicons_file: a tab-delimited text amplicons description file with up to 3 columns, required=True  
  [AMPLICON_NAME, AMPLICON_SEQ, gRNA_SEQ_without_PAM(Optional)]

--threads: number of threads to run sentieon minimap2 & driver module, default=1 
  
--anno: annotate variants with ANNOVAR [1] or not run ANNOVAR [2], required = False  

--assembly: assembly version, hg19,hg38 ...  

--db: ANNOVAR database path, required = False    

--ClinVar: only organism homo sapiens experiment type sequencing data support variant annotations from ClinVar[1], default = 0   

--cleavage_offset: Center of quantification window to use within respect to the 3-end of the provided sgRNA sequence, default = -3  

--window_size: defines the size (in bp) of the quantification window extending from the position specified by the cleavage_offset parameter in relation to the provided guide RNA sequence, 0 means whole amplicon analysis, default = 0  

--o: output path, default='.', required = False  

--ignore_substitutions: enable substitutions evaluation[1], default = 0  

--min_tumor_allele_frac: the minimum allelic fraction in treated sample, default = 0.005  

--min_num_of_reads: the minimum number of reads (per locus site) to evaluate, default = 500  

--max_fisher_pv_active: the maximum pvalue of the statistical difference between treated and untreated sample, default = 0.05  

### 2. Whole genome sequencing (WGS) data analysis
python scripts/CRISPRdetectorWGS.py

--e1: treatment group fq1 path, required = True  

--e2: treatment group fq2 path, required = False  

--c1: control group fq2 path, required = False  

--c2: control group fq2 path, required = False  

--sample: sample name & output dir name    

--o: output path, default='.', required = False  

--threads: number of threads to run sentieon minimap2 & driver module, default = 1   

--bed: bed format file input to call variants of interest region, required = False   

--assembly: path to assembly in fasta format : hg38.fa mm9.fa ... required = True  



### 3. Vector sequence insertion locations detection 
python scripts/CRISPRdetectorVEC.py

--e1: treatment group fq1 path, required = True  

--e2: treatment group fq2 path, required = False  

--c1: control group fq2 path, required = False  

--c2: control group fq2 path, required = False  

--sample: sample name & output dir name  

--o: output path, default='.', required = False  

--threads: number of threads to run sentieon minimap2 & driver module, default = 1   

--bed: bed format file input to call variants of interest region, required = False   

--vector : path to vector genome in fasta format, required = True  

--assembly: path to assembly in fasta format : hg38.fa mm9.fa ... required = True  
