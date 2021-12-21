#$ -q long-sl7
#$ -N tapestri_picard
#$ -e /users/lvelten/project/Methylome/data/stderr.txt #The file that stderr will get written to
#$ -o /users/lvelten/project/Methylome/data/stdout.txt #The file that stdout will get written to
#$ -l virtual_free=10G #Amount of memory
#$ -l h_rt=12:00:00 #Wall time - after this amount of time your job will die. The higher the value the longer you might have to wait in the queue...

echo "###################################################################################################"
echo "$(date): Running PicardTools"
source /software/as/el7.2/anaconda3/bin/activate /users/lvelten/mscherer/conda/envs/tapestri
in_file=/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/bam/Sample8_70_percent_good_performance_aligned.bam
out_file=/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/bam/Sample8_70_percent_good_performance_aligned_fixed.bam
tmp_folder=/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/tmp/
picard FixMateInformation -Xmx4g -Djava.io.tmpdir=$tmp_folder I=$in_file O=$out_file SORT_ORDER=coordinate

