
#!/bin/bash

# remember to chmod +x

# for pair-end reads

salmon quant -i ~/mks/data_and_results/step9_salmon_quant/salmon_index_sa_our/salmon_index_sa \
-l A \
-1 <(zcat SRR9016157.1.fastq.gz) \
-2 <(zcat SRR9016157.2.fastq.gz) \
-p 12 \
-o ~/mks/ouab_theme/step_2_salmon/quants/SRR9016157 \
--numBootstraps 30 \
--seqBias \
--gcBias \
--validateMappings

# for single-end data 

for data in *;
do
	if [[ ${data} == *.fq.gz ]]
	then
	echo "Processing sample ${data}"
	salmon quant -i ~/mks/data_and_results/step9_salmon_quant/salmon_index_sa_our/salmon_index_sa \
             -l A \
             -r <(zcat ${data}) \
             -p 12 -o ~/mks/data_and_results/step9_salmon_quant/quants/$(echo ${data} | cut -c 1-11) \
             --numBootstraps 30 \
             --seqBias \
             --gcBias \
             --validateMappings
	else
	echo "File is not reads data"
	fi
done

