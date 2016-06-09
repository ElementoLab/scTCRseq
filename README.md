# scTCRseq
## Introduction

This project is an implementation of a pipeline for Single-cell RNAseq package for recovering TCR data in python

[Github Project](https://github.com/ElementoLab/scTCRseq)

## Configuration and Dependencies
The pipeline needs for the following programs to be installed and the paths :

###SEQTK:
https://github.com/lh3/seqtk

###Blastall:
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

###GapFiller:
http://www.baseclear.com/genomics/bioinformatics/basetools/gapfiller

###Vidjil:
https://github.com/vidjil/vidjil

###And their accompanying paths need to be changed in the script cmd_line_sctcrseq.py:

seqTkDir="/path/to/seqtk/"

blastallDir="/path/to/blastall/"

gapFillerDir="/path/to/GapFiller_v1-10_linux-x86_64/"

vidjildir="/path/to/vidjil/"

lengthScript="/path/to/calc.median.read.length.pl"

###Also the user can select their chosen TCR alpha and beta V and C reference databases (we recommend downloading from imgt.org) and enter their locations:

location for FASTA BLAST reference sequences downloadable from imgt.org - (needs to be changed manually)


humanTRAVblast="/path/to/TRAV.human.fa"

humanTRBVblast="/path/to/TRBV.human.fa"

humanTRACblast="/path/to/TRAC.human.fa"

humanTRBCblast="/path/to/TRBC.human.fa"

mouseTRAVblast="/path/to/TRAV.mouse.fa"

mouseTRBVblast="/path/to/TRBV.mouse.fa"

mouseTRACblast="/path/to/TRAC.mouse.fa"

mouseTRBCblast="/path/to/TRBC.mouse.fa"



###location for Vidjil BLAST reference sequences in vidjil program -  (needs to be changed manually)

humanVidjilRef="/path/to/tr_germline/human"

mouseVidjilRef="/path/to/tr_germline/mouse"



## Example Command Line

We recommend running the pipeline on paired end fluidigm single cell RNA seq data.

the usage is as follows:

python cmd_line_sctcrseq.py --fastq1 FASTQ1 --fastq2 FASTQ2 --species human/mouse --outdir OUTPUT DIRECTORY --label OUTPUT LABEL

also running

python cmd_line_sctcrseq.py 

Will give command line options


to test the program we recommend either using your own single cell RNA sequencing data or downloading data such as at 

http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1104129



 


