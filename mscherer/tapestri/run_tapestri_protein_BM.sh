#$ -q long-sl7
#$ -N tapestri_protein_Sample12
#$ -e /users/lvelten/project/Methylome/analysis/missionbio/tapestri/protein/Sample12/stderr.txt #The file that stderr will get written to
#$ -o /users/lvelten/project/Methylome/analysis/missionbio/tapestri/protein/Sample12/stdout.txt #The file that stdout will get written to
#$ -l virtual_free=30G #Amount of memory
#$ -l h_rt=96:00:00 #Wall time - after this amount of time your job will die. The higher the value the longer you might have to wait in the queue...
#$ -pe smp 4 

/users/lvelten/project/Methylome/src/scTAM-seq-scripts/mscherer/tapestri/tapestri_pipeline_protein_BCells.sh -p /users/lvelten/project/Methylome/data/BM/Sample12/protein/ -s Sample12_protein -n Sample12 -o /users/lvelten/project/Methylome/analysis/missionbio/tapestri/protein/Sample12/ -c 4
