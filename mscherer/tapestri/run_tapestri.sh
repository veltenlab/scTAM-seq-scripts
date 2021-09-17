#$ -q long-sl7
#$ -N tapestri_Sample5_80_1_percent
#$ -e /users/lvelten/project/Methylome/analysis/missionbio/downsampling/Sample5/stderr.txt #The file that stderr will get written to
#$ -o /users/lvelten/project/Methylome/analysis/missionbio/downsampling/Sample5/stdout.txt #The file that stdout will get written to
#$ -l virtual_free=100G #Amount of memory
#$ -l h_rt=96:00:00 #Wall time - after this amount of time your job will die. The higher the value the longer you might have to wait in the queue...
#$ -pe smp 16 

/users/lvelten/project/Methylome/src/scTAM-seq-scripts/mscherer/tapestri/tapestri_pipeline.sh -p /users/lvelten/project/Methylome/data/raw/ -s Sample5_downsampled -n Sample5_70_percent -o /users/lvelten/project/Methylome/analysis/missionbio/downsampling/Sample5_70_percent/ -c 16
