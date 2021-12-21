#$ -q long-sl7
#$ -N tapestri_bwa
#$ -e /users/lvelten/project/Methylome/data/stderr.txt #The file that stderr will get written to
#$ -o /users/lvelten/project/Methylome/data/stdout.txt #The file that stdout will get written to
#$ -l virtual_free=10G #Amount of memory
#$ -l h_rt=48:00:00 #Wall time - after this amount of time your job will die. The higher the value the longer you might have to wait in the queue...
#$ -pe smp 16

echo "###################################################################################################"
echo "$(date): Alignment with BWA"
in_first=/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/barcode/Sample8_70_percent_good_performance_read1.with.barcodes.fastq.gz
in_second=/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/barcode/Sample8_70_percent_good_performance_read2.with.barcodes.fastq.gz
out_file=/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance//bam/Sample8_70_percent_good_performance_aligned.bam
sam_file=/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/tmp//tmp.sam
tmp_bam=/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/tmp//tmp.bam
cores=4
ref_genome="/users/lvelten/project/Methylome/references/MissionBio/hg19/ucsc_hg19.fa"
sam_header=/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/tmp/Sample8_70_percent_good_performance.sam.header
bwa mem -C -M -t $cores -H $sam_header $ref_genome $in_first $in_second > $sam_file
samtools view -@ $cores -S -b $sam_file > $tmp_bam
samtools view -@ $cores -bh -q 30 -F 4 -F 8 -F 0X0100 $tmp_bam > $out_file
