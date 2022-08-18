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

else
    mkdir $output
fi

echo "###################################################################################################"
echo "$(date): Start MissionBio Tapestri Pipeline"
# Variables
#adapter_sequence_a="CTGTCTCTTATA"
#adapter_sequence_g="TATAAGAGACAG"
#To be downloaded from: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
ref_genome="ucsc_hg19.fa"
panel_bed="panel.bed"
ampli_file="Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv"
cellfinder_cutoff=0.7
tmp_folder=${output}/tmp
mkdir $tmp_folder

# Load conda env
conda activate $TAPESTRI

# PART 1:

# Barcodecode (i.e. extracting the barcodes from the fastqs
echo "###################################################################################################"
echo "$(date): Trimming and Barcode Extraction"
barcode_folder=${output}/barcode
mkdir $barcode_folder
in_first=${path}/${sample}_read1.fastq.gz
in_second=${path}/${sample}_read2.fastq.gz
out_first=${barcode_folder}/${name}_read1.with.barcodes.fastq.gz
out_second=${barcode_folder}/${name}_read2.with.barcodes.fastq.gz
barcode_output=${barcode_folder}/${name}_raw.barcodes.txt
barcode_metrics=${barcode_folder}/${name}_01_barcode.metrics.json
sam_header=${tmp_folder}/${name}.sam.header
tapestri barcode extract --input-fwd $in_first --input-rev $in_second --output-fwd $out_first --output-rev $out_second --output-barcodes $barcode_output --output-metrics $barcode_metrics --chemistry dna.V2 --cut-adapters --decode --cores $cores --barcode-suffix -1
sed -n 's/^\([^\t]*\).*/@RG\tID:\1\tSM:\1/p' $barcode_output > $sam_header

## BWA alignment
echo "###################################################################################################"
echo "$(date): Alignment with BWA"
in_first=$out_first
in_second=$out_second
out_file=${output}/bam/
mkdir $out_file
out_file=${out_file}/${name}_aligned.bam
sam_file=${tmp_folder}/tmp.sam
tmp_bam=${tmp_folder}/tmp.bam
bwa mem -C -M -t $cores -H $sam_header $ref_genome $in_first $in_second > $sam_file
samtools view -@ $cores -S -b $sam_file > $tmp_bam
samtools view -@ $cores -bh -q 30 -F 4 -F 8 -F 0X0100 $tmp_bam > $out_file
#rm -rf $sam_file
#rm -rf $tmp_bam

## Running picard
echo "###################################################################################################"
echo "$(date): Running PicardTools"
in_file=$out_file
out_file=${output}/bam/${name}_aligned_fixed.bam
picard FixMateInformation -Xmx4g -Djava.io.tmpdir=$tmp_folder I=$in_file O=$out_file SORT_ORDER=coordinate

## Computing the barcode distribution
echo "###################################################################################################"
echo "$(date): Computing the barcode distribution"
in_file=$out_file
count_file=${output}/${name}.mapped.target.count.txt
out_file=${output}/${name}.cellfinder.barcode.distribution.txt
minus_5=${tmp_folder}/${name}_panel_minus_5.bed
awk '{print $1 "\t" ($2 + 5) "\t" ($3 - 5) "\t" $4}' $panel_bed > $minus_5
samtools view -@ $cores $in_file | grep -oP 'RG:Z:\K.+?\b' | sort | uniq -c | sort -n -r | awk '{$1=$1}1' > $count_file
perl $TAPESTRI/bin/resources/Perl/distribution_barcodes_amplicon.pl --amplicon-input $minus_5 --bam-input $in_file --output-r1 $out_file --output-r2 ${output}/${name}.r2.barcode.distribution.txt

# Running the R-script resembling CellFinder (our version)
echo "###################################################################################################"
echo "$(date): Running our version of CellFinder"
conda activate $RNBEADS
in_file=$out_file
out_file=${output}/tsv/
mkdir $out_file
out_file=${out_file}/${name}.barcode.cell.distribution.tsv
Rscript cellfinder_BCells_selected.R -f $in_file -a $ampli_file -c $cellfinder_cutoff -o $out_file

# Calculate reads per cell    
echo "###################################################################################################"
echo "$(date): Extracting the reads mapped to cells"
conda activate $TAPESTRI
in_file=$out_file
out_file=${tmp_folder}/cells.txt
awk '{if (NR > 1) print $1}' $in_file > $out_file
in_file=$out_file
out_file=${output}/bam/${name}.cells.bam
samtools view -@ $cores -b -R $in_file ${output}/bam/${name}_aligned_fixed.bam > $out_file
samtools index $out_file
python3 -m missionbio.dna.calculate_reads_mapped_to_cells --tsv ${output}/tsv/${name}.barcode.cell.distribution.tsv > ${output}/tsv/read_to_cells.txt

# Run DoubletDetection
echo "###################################################################################################"
echo "$(date): Running DoubletDetection"
in_doublet=${output}/tsv/${name}.barcode.cell.distribution.tsv
out_doublet=${output}/tsv/doublet_scores_DoubletDetection.csv
python3 DoubletDetection.py --input $in_doublet --output $out_doublet

# DONE!
echo "###################################################################################################"
echo "$(date): Finished PART I"


# DONE!
echo "###################################################################################################"
echo "$(date): Finished PART I"

echo "###################################################################################################"
echo "$(date): Finished the MissionBio pipeline"
