#!/bin/bash
# Run this script from the yeast-rnaseq-pipeline project root folder
set -e
mkdir -p data/raw data/processed results

ACCESSIONS="ERR458493 ERR458494 ERR458495 ERR458500 ERR458501 ERR458502"
REF=data/raw/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
INDEX=data/processed/Saccharomyces.index

if [ ! -f $REF ]; then
   wget https://ftp.ensembl.org/pub/current_fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -P data/raw
fi

#index the reference file

if [ ! -d $INDEX ]; then
   salmon index -t $REF -i $INDEX
fi


#loop through the RNA seq And Quant
for accession in $ACCESSIONS
do
   if [ ! -f data/raw/${accession}.fastq.gz ]; then
        wget --tries=3 --waitretry=30 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/${accession}/${accession}.fastq.gz -P data/raw
   fi

   if [ ! -f results/${accession}_quant/quant.sf ]; then
        salmon quant -i $INDEX -r data/raw/${accession}.fastq.gz -o results/${accession}_quant -l A
   fi
done 

#quatmerge and for counting matrix
salmon quantmerge --quants results/*_quant -c numreads -o data/processed/combined_matrix.tsv

python script/yeast-rnaseq-pipeline.py