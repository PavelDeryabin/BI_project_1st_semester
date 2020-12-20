# Installation of the salmon via conda

https://anaconda.org/bioconda/salmon

# operations for the index generarion with specifik k-mer size equal to 21 

1. downloading last release of the complete h gemone from /pub/databases/gencode/Gencode_human/release_33/

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh38.primary_assembly.genome.fa.gz
    
2. Preparing metadata

grep "^>" <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt 

sed -i.bak -e 's/>//g' decoys.txt

sed -i -e 's/>//g' decoys.txt
    
3. downloading coding and non-coding transcriptome 

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.lncRNA_transcripts.fa.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.transcripts.fa.gz
    
4. building an concatenated transcriptome and genome reference file for index 

cat gencode.v33.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz
    
5. building final gentrome at the different k-mers

salmon index -t gentrome.fa.gz --decoys decoys.txt -p 12 -i salmon_index_sa --gencode -k 21 


# script6.sh describes operations nessessary for transcripts abundances estimation 
