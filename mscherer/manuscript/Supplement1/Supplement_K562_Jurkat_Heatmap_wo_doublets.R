############## Supplement_K562_Jurkat_Heatmap.R ##############
#' This script is used to plot a binary heatmap of the read counts. A read cutoff has to be specified (default 10)
#' to differentiate between presence and absence of read counts, i.e., DNA methylation. Additionally, the read count
#' file and the assignment of cells to clusters has to be specified.

library(viridis)
library(pheatmap)
cut.off <- 0
sample <- 'Sample5_80_percent'
plot.path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Supplement/'
dat.s3 <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'))
to.plot <- apply(dat.s3,2,function(x){
  ifelse(x>cut.off,1,0)
})
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
doublet.score <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores_DoubletDetection.csv'), row.names=2)
ampli.info <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/amplicons/cell_lines/R/amplicon_info_with_RnBeads.csv')
row.names(ampli.info) <- ampli.info$AmpID_design_1898
colnames(to.plot) <- colnames(dat.s3)
row.names(to.plot) <- row.names(dat.s3)
clust.file <- clust.file[row.names(to.plot),]
row.names(clust.file) <- row.names(to.plot)
clust.file$Doublet <- ifelse(doublet.score[row.names(clust.file), 'DoubletDetectionLabel']==0, 'Singlet', 'Doublet')
to.plot <- to.plot[which(clust.file$Doublet=='Singlet'),]
ampli.info$Type[grepl('Aci',ampli.info$Type)] <- 'Aci'
ampli.info$Type <- as.factor(ampli.info$Type)
to.plot <- to.plot[,row.names(ampli.info)]
to.plot <- to.plot[,order(ampli.info$Type)]
ampli.info <- ampli.info[order(ampli.info$Type),]
#to.plot <- to.plot[!(clust.file$CellType%in%'Mixed'),]
for(cat in unique(ampli.info$Type)){
  sel.cols <- ampli.info$Type%in%cat
  sel.means <- apply(to.plot[,sel.cols],2,mean,na.rm=TRUE)
  to.plot[,sel.cols] <- to.plot[,sel.cols][,order(sel.means)]
}
png(file.path(plot.path, 'Supplement_K562_Jurkat_Heatmap_wo_Doublets.png'))
pheatmap(to.plot,
         annotation_row = clust.file[,c('CellType'), drop=FALSE],
         color=rev(inferno(50)),
         annotation_col = ampli.info[,c('Type'),drop=FALSE],
         show_rownames=FALSE,show_colnames=FALSE,
         clustering_distance_cols = "binary",
         clustering_method = 'ward.D2',
         cluster_cols = FALSE,
         annotation_colors=list(CellType=c('Jurkat'='#ff9289', 'K562'='#82b7ff', 'Mixed'='grey'),
                                Type=c('Aci'='grey','Hha Jurkat high'='#b3645eff','Hha K562 high'='#597daeff'),
                                Doublet=c('Singlet'='white','Doublet'='black')),
         legend=FALSE,
         annotation_legend=FALSE)
dev.off()
