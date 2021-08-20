#$ -q short-sl7
#$ -N tapestri_Sample3_80
#$ -e /users/lvelten/project/Methylome/analysis/missionbio/tapestri/protein/Sample3_80_percent/stderr.txt #The file that stderr will get written to
#$ -o /users/lvelten/project/Methylome/analysis/missionbio/tapestri/protein/Sample3_80_percent/stdout.txt #The file that stdout will get written to
#$ -l virtual_free=100G #Amount of memory
#$ -l h_rt=06:00:00 #Wall time - after this amount of time your job will die. The higher the value the longer you might have to wait in the queue...
#$ -pe smp 16 

/users/lvelten/project/Methylome/src/scTAM-seq-scripts/mscherer/tapestri/tapestri_pipeline_protein.sh -p /users/rbeekman/abianchi/Mission_Bio/Cell_lines_data_2/BEEKMANREN_05 -s T5PROT_lib_07978AAC_CAGAGAGG-ACGAGAAC -n Sample5_80_percent -o /users/lvelten/project/Methylome/analysis/missionbio/tapestri/protein/Sample5_80_percent/protein/ -c 16
