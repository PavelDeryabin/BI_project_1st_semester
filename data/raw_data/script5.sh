
#!/bin/bash

for data in *;
do
	if [[ ${data} == *.fastq.gz ]]
	then
	echo "Processing sample ${data}"
	salmon quant -i ~/mks/data_and_results/step9_salmon_quant/salmon_index_sa_our/salmon_index_sa \
             -l A \
             -r <(zcat ${data}) \
             -p 12 \
	     -o ~/mks/ouab_theme/step_2_salmon/quants/$(echo ${data} | cut -c 1-10) \
             --numBootstraps 30 \
             --seqBias \
             --gcBias \
             --validateMappings
	else
	echo "File is not reads data"
	fi
done


