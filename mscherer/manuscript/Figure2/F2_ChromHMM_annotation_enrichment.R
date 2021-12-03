######################## F2_ChromHMM_annotation_enrichment.R ######################## 
#' This file generates enrichment analyses for the chromatin states identified in 10.1016/j.ccell.2016.09.014
#' for the variable CpGs that we identified for each of the clusters.

library(RnBeads)
all.cpgs <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
all.cpgs <- subset(all.cpgs, subset = Type.of.amplicon=='CpG.B.cell.diff', select='background.cpgs')
background.cpgs <- all.gr[all.cpgs[, 'background.cpgs']]
cell_metadata <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/rowinfo_reclustering_doublet.csv',
                          row.names=1)
file_map <- c('Cluster1'='/users/mscherer/cluster/project/Methylome/infos/BCells/chromatin_states/NBCB_12_segments_intersect.bed',
              'Cluster2a'='/users/mscherer/cluster//project/Methylome/infos/BCells/chromatin_states/ncsMBC_12_segments.bed',
              'Cluster2b'='/users/mscherer/cluster//project/Methylome/infos/BCells/chromatin_states/csMBC_12_segments_intersect.bed',
              'Cluster2c'='/users/mscherer/cluster//project/Methylome/infos/BCells/chromatin_states/csMBC_12_segments_intersect.bed')
res <- list()
for(clust in names(file_map)){
  cluster.cpgs <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/amplicon_correlation/variable_amplicons_Cluster1.csv')
  variable.cpgs <- background.cpgs[cluster.cpgs$CpGID]
  anno <- read.table(file_map[clust])
  anno <- makeGRangesFromDataFrame(anno,
                                   seqnames.field = 'V1',
                                   start.field = 'V2',
                                   end.field = 'V3',
                                   keep.extra.columns=TRUE)
  part_res <- list()
  for(state in unique(values(anno)[, 1])){
    sel_anno <- anno[values(anno)[, 1]%in%state]
    tps <- length(findOverlaps(variable.cpgs, sel_anno))
    fps <- length(variable.cpgs)-tps
    fns <- length(findOverlaps(background.cpgs, sel_anno))
    tns <- length(background.cpgs)-fns
    gr <- fisher.test(matrix(c(tps,fns,fps,tns),2,2),alternative="greater")$p.value
    or <- (tps/fps)/(fns/tns)
    part_res[[state]] <- c(OR=or, PValue=gr)
  }
  res[[clust]] <- part_res
}