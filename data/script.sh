#!/bin/bash

# Example for SRR5931953 (A549 - ctr - GSE102639)

# Download SRA files from NCBI GEO with SRA-toolkit

prefetch SRR5931953 

# Convert to .fastq

fasterq-dump --split-files SRR5931953.sra
rm SRR5931953.sra # removing .sra files









































