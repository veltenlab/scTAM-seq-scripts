#!/bin/bash

while getopts p:s:n:o:c: flag
do
    case "${flag}" in
        p) path=${OPTARG};;
        s) sample=${OPTARG};;
        n) name=${OPTARG};;
        o) output=${OPTARG};;
        c) cores=${OPTARG};;
    esac
done

if [ -f $output ]
then
   echo "Output directory must not exist"; exit $ERRCODE;
fi
#else
#    mkdir $output
#fi

echo "###################################################################################################"
echo "$(date): Started the MissionBio Protein pipeline"
# Variables
#adapter_sequence_a="CTGTCTCTTATA"
#adapter_sequence_g="TATAAGAGACAG"
ref_genome="/users/lvelten/project/Methylome/references/MissionBio/hg19/ucsc_hg19.fa"
adapter_3="/users/lvelten/project/Methylome/data/MissionBio/protein_pipeline/ab_adapters_3.fasta"
adapter_5="/users/lvelten/project/Methylome/data/MissionBio/protein_pipeline/ab_adapters_5.fasta"
ab_barcodes="/users/lvelten/project/Methylome/data/MissionBio/protein_pipeline/PG-135.csv"

# Load conda env
source /software/as/el7.2/anaconda3/bin/activate /users/lvelten/mscherer/conda/envs/tapestri

# Run the pipeline
echo "###################################################################################################"
echo "$(date): Starting the protein pipeline command"
in_first=${path}/${sample}_R1_001.fastq.gz
in_second=${path}/${sample}_R2_001.fastq.gz
tapestri protein run --ab-adapters-3 $adapter_3 --ab-adapters-5 $adapter_5 --ab-barcodes $ab_barcodes --output-folder $output --output-prefix $name --n-cores $cores --r1 $in_first --r2 $in_second

# Cleanup
echo "###################################################################################################"
echo "$(date): Cleanup the mess"
fastq=${output}/results/trimmed-fastq/
reads=${output}/results/tsv/${name}-ab-reads.tsv.gz
rm -rf $fastq
rm -rf $reads

echo "###################################################################################################"
echo "$(date): Finished the MissionBio Protein pipeline"
