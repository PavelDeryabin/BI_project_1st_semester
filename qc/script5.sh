
#!/bin/bash

for data in *;
do
	if [[ ${data} == *.fastq.gz ]]
	then
	echo "Processing sample ${data}"
	trimFilter --ifq ${data} \
           --length 50 \
           --output trimmed/$(echo ${data} | cut -c 1-6) \
           --gzip n \
           --trimQ ENDS \
           --minL 25 \
           --minQ 30 \
           --zeroQ 33 \
           --percent 5 \
           --trimN ENDS
	cd trimmed
	gzip $(echo ${data} | cut -c 1-6)_good.fq
	cd ..
	else
	echo "File is not reads data"
	fi
done

cd trimmed
rm *_NNNN.fq
rm *_cont.fq
rm *_lowq.fq

