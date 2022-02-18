#$ -q long-sl7
#$ -N JOBNAME
#$ -e OUTFOLDER/stderr.txt #The file that stderr will get written to
#$ -o OUTFOLDER/stdout.txt #The file that stdout will get written to
#$ -l virtual_free=40G #Amount of memory
#$ -l h_rt=96:00:00 #Wall time - after this amount of time your job will die. The higher the value the longer you might have to wait in the queue...
#$ -pe smp 16 

#INPATH=folder which stores the raw fastq data
#INNAME=name of the sample proceeding the _read1.fastq.gz and _read2.fastq.gz
#OUTNAME=name of the output files
#OUTFOLDER=the output folder
tapestri_pipeline_BCells.sh -p $INPATH -s $INNAME -n $OUTNAME -o $OUTFOLDER -c 16
