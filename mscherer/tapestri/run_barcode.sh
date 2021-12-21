#$ -q short-sl7
#$ -N tapestri_barcode
#$ -e /users/lvelten/project/Methylome/data/stderr.txt #The file that stderr will get written to
#$ -o /users/lvelten/project/Methylome/data/stdout.txt #The file that stdout will get written to
#$ -l virtual_free=10G #Amount of memory
#$ -l h_rt=05:55:00 #Wall time - after this amount of time your job will die. The higher the value the longer you might have to wait in the queue...

## Computing the barcode distribution
echo "###################################################################################################"
echo "$(date): Computing the barcode distribution"
source /software/as/el7.2/anaconda3/bin/activate /users/lvelten/mscherer/conda/envs/tapestri
in_file=/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/bam/Sample8_70_percent_good_performance_aligned_fixed.bam
count_file=/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/Sample8_70_percent_good_performance.mapped.target.count.txt
out_file=/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/Sample8_70_percent_good_performance.cellfinder.barcode.distribution.txt
minus_5=/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/tmp/Sample8_70_percent_good_performance_panel_minus_5.bed
panel_bed="/users/lvelten/project/Methylome/infos/BCells/panel.bed"
cores=1
awk '{print $1 "\t" ($2 + 5) "\t" ($3 - 5) "\t" $4}' $panel_bed > $minus_5
#samtools view -@ $cores $in_file | grep -oP 'RG:Z:\K.+?\b' | sort | uniq -c | sort -n -r | awk '{$1=$1}1' > $count_file
perl /users/lvelten/mscherer/conda/envs/tapestri/bin/resources/Perl/distribution_barcodes_amplicon.pl --amplicon-input $minus_5 --bam-input $in_file --output-r1 $out_file --output-r2 /users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/Sample8_70_percent_good_performance.r2.barcode.distribution.txt

