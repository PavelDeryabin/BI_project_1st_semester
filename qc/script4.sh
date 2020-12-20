#!/bin/bash

for data in *.fastq.gz;
do
echo "Processing sample ${data}"
~/Software/BBMap_38.75/bbmap/filterbytile.sh in=${data} out=filtered/${data} 
done

