#!/bin/bash

# for pair-end data

~/Software/BBMap_38.75/bbmap/bbduk.sh in1=data_1.fastq.gz in2=data_2.fastq.gz out1=res_1.fastq.gz out2=res_2.fastq.gz ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo ftl=10

# for single-end data

~/Software/BBMap_38.75/bbmap/bbduk.sh in=data.fastq.gz out=res.fastq.gz ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo ftl=10
 
