#$ -q long-sl7
#$ -N tapestri_protein_Sample12
#$ -e $OUTPUT/stderr.txt #The file that stderr will get written to
#$ -o $OUTPUT/stdout.txt #The file that stdout will get written to
#$ -l virtual_free=30G #Amount of memory
#$ -l h_rt=96:00:00 #Wall time - after this amount of time your job will die. The higher the value the longer you might have to wait in the queue...
#$ -pe smp 4 

tapestri_pipeline_protein_BCells.sh -p $INPUT -s Sample12 -n Sample12 -o $OUTPUT -c 4
