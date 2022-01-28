#$ -q short-sl7
#$ -N tapestri_cellfinder
#$ -e /users/lvelten/project/Methylome/data/stderr.txt #The file that stderr will get written to
#$ -o /users/lvelten/project/Methylome/data/stdout.txt #The file that stdout will get written to
#$ -l virtual_free=50G #Amount of memory
#$ -l h_rt=00:30:00 #Wall time - after this amount of time your job will die. The higher the value the longer you might have to wait in the queue...

echo "###################################################################################################"
echo "$(date): Running our version of CellFinder"
source /software/as/el7.2/anaconda3/bin/activate /users/lvelten/mscherer/conda/envs/rnbeads
in_file=/users/lvelten/project/Methylome/analysis/missionbio/BM/Sample12_70_percent_good_performance/Sample12_70_percent_good_performance.cellfinder.barcode.distribution.txt
out_file=/users/lvelten/project/Methylome/analysis/missionbio/BM/Sample12_70_percent/
mkdir $out_file
out_file=/users/lvelten/project/Methylome/analysis/missionbio/BM/Sample12_70_percent/tsv/
mkdir $out_file
out_file=${out_file}/Sample12_70_percent.barcode.cell.distribution.tsv
ampli_file="/users/lvelten/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv"
cellfinder_cutoff=0.7
Rscript /users/lvelten/project/Methylome/src/scTAM-seq-scripts/mscherer/tapestri/cellfinder_BCells.R -f $in_file -a $ampli_file -c $cellfinder_cutoff -o $out_file

