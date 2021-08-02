#$ -q long-sl7
#$ -N tapestri_Sample3_80
#$ -e /users/lvelten/project/Methylome/analysis/missionbio/tapestri/Sample3_80_percent/stderr.txt #The file that stderr will get written to
#$ -o /users/lvelten/project/Methylome/analysis/missionbio/tapestri/Sample3_80_percent/stdout.txt #The file that stdout will get written to
#$ -l virtual_free=100G #Amount of memory
#$ -l h_rt=96:00:00 #Wall time - after this amount of time your job will die. The higher the value the longer you might have to wait in the queue...
#$ -pe smp 16 

/users/lvelten/project/Methylome/src/MissionBio/tapestri_pipeline.sh -p /users/rbeekman/johernandez/backup-sequencing/210201_NS500645_0245_AHK5N7BGXH/ -s sDcKbn210201ae -n Sample3_80_percent -o /users/lvelten/project/Methylome/analysis/missionbio/tapestri/Sample3_80_percent/ -c 16
