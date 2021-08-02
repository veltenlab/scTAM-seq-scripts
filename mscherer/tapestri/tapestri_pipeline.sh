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
ref_genome="/users/lvelten/project/Methylome/references/MissionBio/hg19/ucsc_hg19.fa"
panel_bed="/users/lvelten/project/Methylome/infos/cell_lines/TapestriDesigner/1898.bed"
ampli_file="/users/lvelten/project/Methylome/infos/amplicon_info/200804_Design_Mission_Bio_sites_Jurkat_K562.tsv"
cellfinder_cutoff=0.8
tmp_folder=${output}/tmp
mkdir $tmp_folder

# Load conda env
source /software/as/el7.2/anaconda3/bin/activate /users/lvelten/mscherer/conda/envs/tapestri

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
rm -rf $sam_file
rm -rf $tmp_bam

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
perl /users/lvelten/mscherer/conda/envs/tapestri/bin/resources/Perl/distribution_barcodes_amplicon.pl --amplicon-input $minus_5 --bam-input $in_file --output-r1 $out_file --output-r2 ${output}/${name}.r2.barcode.distribution.txt

# Running the R-script resembling CellFinder (our version)
echo "###################################################################################################"
echo "$(date): Running our version of CellFinder"
source /software/as/el7.2/anaconda3/bin/activate /users/lvelten/mscherer/conda/envs/rnbeads
in_file=$out_file
out_file=${output}/tsv/
mkdir $out_file
out_file=${out_file}/${name}.barcode.cell.distribution.tsv
Rscript /users/lvelten/project/Methylome/src/MissionBio/cellfinder.R -f $in_file -a $ampli_file -c $cellfinder_cutoff -o $out_file

# Calculate reads per cell    
echo "###################################################################################################"
echo "$(date): Extracting the reads mapped to cells"
source /software/as/el7.2/anaconda3/bin/activate /users/lvelten/mscherer/conda/envs/tapestri
in_file=$out_file
out_file=${tmp_folder}/cells.txt
awk '{if (NR > 1) print $1}' $in_file > $out_file
in_file=$out_file
out_file=${output}/bam/${name}.cells.bam
samtools view -@ $cores -b -R $in_file ${output}/bam/${name}_aligned_fixed.bam > $out_file
samtools index $out_file
python3 -m missionbio.dna.calculate_reads_mapped_to_cells --tsv ${output}/tsv/${name}.barcode.cell.distribution.tsv > ${output}/tsv/read_to_cells.txt

# DONE!
echo "###################################################################################################"
echo "$(date): Finished PART I"

# Genotype calling:

# Split the cells.bam file by barcode
echo "###################################################################################################"
echo "$(date): Split bamfile by barcode"
in_file=$out_file
out_folder=${tmp_folder}/bam/
mkdir $out_folder
python3 /users/lvelten/mscherer/conda/envs/tapestri/lib/python3.6/site-packages/missionbio/dna/split_bam_by_barcode.py -i $in_file -o $out_folder -n $cores

# Running GATK HaplotypeCaller
echo "###################################################################################################"
echo "$(date): Start the HaplotypeCaller"
in_folder=$out_folder
out_folder=${tmp_folder}/vcf/
mkdir $out_folder
python3 /users/lvelten/mscherer/conda/envs/tapestri/lib/python3.6/site-packages/missionbio/dna/gatk_haplotype_caller.py --input-dir $in_folder --output-dir $out_folder --config /users/lvelten/project/Methylome/src/MissionBio/config.yaml --prefix gatk

# Combine VCF files using the Genomics import and then combination
echo "###################################################################################################"
echo "$(date): Combine the VCF files"
in_folder=$out_folder
variant_file=${tmp_folder}/variants.txt
db_folder=${tmp_folder}/gendb/
out_file=${output}/vcf/
mkdir $out_file
out_file=${out_file}/${name}.combined.vcf
call_string=""
for f in ${in_folder}/*.vcf.gz; do
    call_string+=" --variant ${f}"
done

gatk_dir='/users/lvelten/mscherer/conda/envs/tapestri/bin/gatk'
gatk GenomicsDBImport $call_string --genomicsdb-workspace-path $db_folder --intervals $panel_bed --merge-input-intervals --batch-size 150
gatk GenotypeGVCFs --java-options -Xmx55g -R $ref_genome --heterozygosity 0.001 --max-alternate-alleles 2 -V gendb://${db_folder} -L $panel_bed --all-sites --include-non-variant-sites --genomicsdb-use-vcf-codec -O $out_file
bgzip $out_file
tabix ${out_file}.gz

# PART II

# Generate the h5 file for further exploration
echo "##########/users/lvelten/project/Methylome/src/MissionBio/metadata.json#########################################################################################"
echo "$(date): Generate the final h5 file"
in_file=${out_file}.gz
count_file=${output}/tsv/${name}.barcode.cell.distribution.tsv
out_file=${output}/${name}.dna.h5
tapestri h5 create dna \
--vcf $in_file \
--read-counts $count_file \
--metadata /users/lvelten/project/Methylome/src/MissionBio/metadata.json \
--output $out_file

#Cleaning up
rm -rf $tmp_folder
rm -rf ${output}/bam/${name}_aligned*
rm -rf ${output}/barcode/

echo "###################################################################################################"
echo "$(date): Finished the MissionBio pipeline"
