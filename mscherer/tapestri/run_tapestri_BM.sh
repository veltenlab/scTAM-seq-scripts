#$ -q long-sl7
#$ -N tapestri_Sample12_70_percent_good_performance
#$ -e /users/lvelten/project/Methylome/analysis/missionbio/BM/Sample12_70_percent_good_performance/stderr.txt #The file that stderr will get written to
#$ -o /users/lvelten/project/Methylome/analysis/missionbio/BM/Sample12_70_percent_good_performance/stdout.txt #The file that stdout will get written to
#$ -l virtual_free=40G #Amount of memory
#$ -l h_rt=96:00:00 #Wall time - after this amount of time your job will die. The higher the value the longer you might have to wait in the queue...
#$ -pe smp 16 

/users/lvelten/project/Methylome/src/scTAM-seq-scripts/mscherer/tapestri/tapestri_pipeline_BCells.sh -p /users/lvelten/project/Methylome/data/BM/Sample12/ -s Sample12 -n Sample12_70_percent_good_performance -o /users/lvelten/project/Methylome/analysis/missionbio/BM/Sample12_70_percent_good_performance/ -c 16
