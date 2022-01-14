######################## F2_ChromHMM_annotation_enrichment.R ######################## 
#' This file generates enrichment analyses for the chromatin states identified in 10.1016/j.ccell.2016.09.014
#' for the variable CpGs that we identified for each of the clusters.

.libPaths(c(.libPaths(), '/users/mscherer/conda/envs/rnbeads/lib/R/library/'))
library(RnBeads)
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/Sample8/differential/'
all.cpgs <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
all.cpgs <- subset(all.cpgs, subset = Type.of.amplicon=='CpG.B.cell.diff', select='background.cpgs')
all.gr <- makeGRangesFromDataFrame(rnb.annotation2data.frame(rnb.get.annotation('probes450')))
background.cpgs <- all.gr[all.cpgs[, 'background.cpgs']]
cell_metadata <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/tsv/rowinfo.csv',
                          row.names=1)
file_map <- c('Cluster1'='/users/mscherer/cluster//project/Methylome/infos/BCells/chromatin_states/NBCB_12_segments_intersect_hg19.bed',
              'Cluster2a'='/users/mscherer/cluster//project/Methylome/infos/BCells/chromatin_states/ncsMBC_12_segments_hg19.bed',
              'Cluster2b'='/users/mscherer/cluster//project/Methylome/infos/BCells/chromatin_states/csMBC_12_segments_intersect_hg19.bed')
res <- list()
filtered.counts <- read.table("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/tsv/Sample8_70_percent_good_performance.barcode.cell.distribution.tsv", row.names = 1, header=T)
#amplicon.info <- read.csv("/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv", row.names = 1)
#filtered.counts <- filtered.counts[row.names(cell_metadata), row.names(amplicon.info)]
amplicon.info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
filtered.counts <- filtered.counts[row.names(cell_metadata), row.names(amplicon.info)[amplicon.info$Type.of.amplicon%in%"CpG.B.cell.diff"]]
more_info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt',
                        row.names = 8)
#background.cpgs <- background.cpgs[more_info[row.names(amplicon.info), 'background.cpgs']]
for(clust1 in names(file_map)){
  sel_cells1 <- cell_metadata$Cluster%in%clust1
  for(clust2 in names(file_map)){
    if(clust1==clust2) next 
    sel_cells2 <- cell_metadata$Cluster%in%clust2
    p.vals <- apply(filtered.counts, 2, function(x){
      wilcox.test(x[sel_cells1], x[sel_cells2])$p.value
    })
    mean.diff <- apply(filtered.counts, 2, function(x){
      c1 <- mean(ifelse(x[sel_cells1]>0, 1, 0))
      c2 <- mean(ifelse(x[sel_cells2]>0, 1, 0))
      c1-c2
    })
    p.vals <- p.adjust(p.vals)
    p.vals[is.na(p.vals)] <- 1
    diff_amplicons <- colnames(filtered.counts)[p.vals<0.01&abs(mean.diff)>0.25]
    diff_cpgs <- more_info[diff_amplicons, 'background.cpgs']
    to_write <- data.frame(Amplicon=diff_amplicons,
                           CpGID=diff_cpgs,
                           PValue=p.vals[p.vals<0.01&abs(mean.diff)>0.25],
                           MeanDiff=mean.diff[p.vals<0.01&abs(mean.diff)>0.25])
    diff_cpgs <- background.cpgs[diff_cpgs]
    write.csv(to_write, file.path(plot_path, paste0('differential_CpGs_', paste0(clust1, 'vs', clust2), '.csv')))
    anno <- read.table(file_map[clust1])
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
    res[[paste0(clust1, 'vs', clust2)]] <- part_res
  }
}
res <- as.data.frame(t(as.data.frame(res)))
res <- res[order(res$PValue), ]
write.csv(res, file.path(plot_path, 'ChromHMM_enrichment_pvals.csv'))
