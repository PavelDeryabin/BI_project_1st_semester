#!/bin/bash

# Download SRA files from NCBI GEO with SRA-toolkit

# cd <...>/BI_project_1st_semester/data/raw_data

prefetch SRR5931953 # A549 - ctr - GSE102639
prefetch SRR5931954 # A549 - ctr - GSE102639 
prefetch SRR5931961 # A549 - sen etoposide - GSE102639 
prefetch SRR5931962 # A549 - sen etoposide - GSE102639 

prefetch SRR8145428 # IMR-90 - ctr - GSE122081
prefetch SRR8145429 # IMR-90 - ctr - GSE122081 
prefetch SRR8145430 # IMR-90 - ctr - GSE122081 
prefetch SRR8145446 # IMR-90 - sen OIS - GSE122081
prefetch SRR8145447 # IMR-90 - sen OIS - GSE122081 
prefetch SRR8145448 # IMR-90 - sen OIS - GSE122081

# find . -maxdepth 1 -name "*.sra" -exec mv {} <...>/BI_project_1st_semester/data/raw_data \;
# rm -r SRR*

fasterq-dump --split-files SRR5931953.sra # convertation to .fastq
rm SRR5931953.sra # removing .sra files
gzip SRR5931953.sra.fastq # each .fastq file

fasterq-dump --split-files SRR5931954.sra
rm SRR5931954.sra
gzip SRR5931954.sra.fastq

fasterq-dump --split-files SRR5931961.sra
rm SRR5931961.sra
gzip SRR5931961.sra.fastq

fasterq-dump --split-files SRR5931962.sra
rm SRR5931962.sra
gzip SRR5931962.sra.fastq

fasterq-dump --split-files SRR5931961.sra
rm SRR5931961.sra
gzip SRR5931961.sra.fastq

fasterq-dump --split-files SRR8145428.sra
rm SRR8145428.sra
gzip SRR8145428.sra.fastq

fasterq-dump --split-files SRR8145429.sra
rm SRR8145429.sra
gzip SRR8145429.sra.fastq

fasterq-dump --split-files SRR8145430.sra
rm SRR8145430.sra
gzip SRR8145430.sra.fastq

fasterq-dump --split-files SRR8145446.sra
rm SRR8145446.sra
gzip SRR8145446.sra.fastq

fasterq-dump --split-files SRR8145447.sra
rm SRR8145447.sra
gzip SRR8145447.sra.fastq

fasterq-dump --split-files SRR8145448.sra
rm SRR8145448.sra
gzip SRR8145448.sra.fastq

fastqc *.gz # raw data QC

#find . -maxdepth 1 -name "*.html" -exec mv {} <...>/BI_project_1st_semester/data/qc/fastqc_of_raw_data \; 
#find . -maxdepth 1 -name "*.zip" -exec mv {} <...>/BI_project_1st_semester/data/qc/fastqc_of_raw_data \;






































