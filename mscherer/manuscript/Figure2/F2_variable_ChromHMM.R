######################## F2_variable_ChromHMM.R ######################## 
#' This file generates enrichment analyses for the chromatin states identified in 10.1016/j.ccell.2016.09.014
#' for the variable CpGs that we identified for each of the clusters.

.libPaths(c(.libPaths(), '/users/mscherer/conda/envs/rnbeads/lib/R/library/'))
library(RnBeads)
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/Sample8/differential/'
all.cpgs <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
all.cpgs <- subset(all.cpgs, subset = Type.of.amplicon=='CpG.B.cell.diff', select='background.cpgs')
all.gr <- makeGRangesFromDataFrame(rnb.annotation2data.frame(rnb.get.annotation('probes450')))
background.cpgs <- all.gr[all.cpgs[, 'background.cpgs']]
cell_metadata <- read.table('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/tsv/rowinfo.tsv',
                          row.names=1)
file_map <- c('Cluster1'='/users/mscherer/cluster//project/Methylome/infos/BCells/chromatin_states/NBCB_12_segments_intersect_hg19.bed',
              'Cluster2a'='/users/mscherer/cluster//project/Methylome/infos/BCells/chromatin_states/ncsMBC_12_segments_hg19.bed',
              'Cluster2b'='/users/mscherer/cluster//project/Methylome/infos/BCells/chromatin_states/csMBC_12_segments_intersect_hg19.bed',
              'Cluster2c'='/users/mscherer/cluster//project/Methylome/infos/BCells/chromatin_states/csMBC_12_segments_intersect_hg19.bed')
res <- list()
filtered.counts <- read.table("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/tsv/Sample8_70_percent_good_performance.barcode.cell.distribution.tsv",
                              row.names = 1, header=T)
cell_metadata <- read.table('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/tsv/rowinfo.tsv',
                            row.names=1)
amplicon.info <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt", row.names = 1)
filtered.counts <- filtered.counts[row.names(cell_metadata), row.names(amplicon.info)]
#amplicon.info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
#filtered.counts <- filtered.counts[row.names(cell_metadata), row.names(amplicon.info)[amplicon.info$Type.of.amplicon%in%"CpG.B.cell.diff"]]
more_info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt',
                        row.names = 8)
background.cpgs <- background.cpgs[more_info[row.names(amplicon.info), 'background.cpgs']]
for(clust in names(file_map)){
  diff_cpgs <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/amplicon_correlation/variable_amplicons_', clust, '.csv'))
  diff_cpgs <- background.cpgs[as.character(diff_cpgs$CpGID)]
  anno <- read.table(file_map[clust])
  anno <- makeGRangesFromDataFrame(anno,
                                   seqnames.field = 'V1',
                                   start.field = 'V2',
                                   end.field = 'V3',
                                   keep.extra.columns=TRUE)
  part_res <- list()
  for(state in unique(values(anno)[, 1])){
    sel_anno <- anno[values(anno)[, 1]%in%state]
    tps <- length(findOverlaps(diff_cpgs, sel_anno))
    fps <- length(diff_cpgs)-tps
    fns <- length(findOverlaps(background.cpgs, sel_anno))
    tns <- length(background.cpgs)-fns
    gr <- fisher.test(matrix(c(tps,fns,fps,tns),2,2),alternative="greater")$p.value
    or <- (tps/fps)/(fns/tns)
    part_res[[state]] <- c(OR=or, PValue=gr)
  }
  res[[clust]] <- part_res
}
res <- as.data.frame(t(as.data.frame(res)))
res <- res[order(res$PValue), ]
