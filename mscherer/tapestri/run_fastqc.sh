#$ -q long-sl7
#$ -N tapestri_fastqc
#$ -e /users/lvelten/project/Methylome/data/stderr.txt #The file that stderr will get written to
#$ -o /users/lvelten/project/Methylome/data/stdout.txt #The file that stdout will get written to
#$ -l virtual_free=10G #Amount of memory
#$ -l h_rt=32:00:00 #Wall time - after this amount of time your job will die. The higher the value the longer you might have to wait in the queue...
#$ -pe smp 1

source /software/as/el7.2/anaconda3/bin/activate /users/lvelten/mscherer/conda/envs/genotools
fastqc /users/lvelten/project/Methylome/data/Sample6/BCells_Sample6_full_read* --outdir /users/lvelten/project/Methylome/data/Sample6/
fastqc /users/lvelten/project/Methylome/data/Sample8/BCells_Sample8_read* --outdir /users/lvelten/project/Methylome/data/Sample8/
fastqc /users/lvelten/project/Methylome/data/Sample8/protein/Sample8_protein_read* --outdir /users/lvelten/project/Methylome/data/Sample8/protein/
