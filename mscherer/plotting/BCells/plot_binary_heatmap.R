############## plot_binary_heatmap.R ##############
#' This script is used to plot a binary heatmap of the read counts. A read cutoff has to be specified (default 10)
#' to differentiate between presence and absence of read counts, i.e., DNA methylation. Additionally, the read count
#' file and the assignment of cells to clusters has to be specified.

library(pheatmap)
library(viridis)
cut.off <- 5
sample <- 'BCells_Sample6_50_percent'
dat.s3 <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'))
to.plot <- apply(dat.s3,2,function(x){
  ifelse(x>cut.off,1,0)
})
ampli.info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.tsv')
colnames(to.plot) <- colnames(dat.s3)
row.names(to.plot) <- row.names(dat.s3)
ampli.info$Type <- as.factor(ampli.info$Type)
to.plot <- to.plot[,row.names(ampli.info)]
to.plot <- to.plot[,order(ampli.info$Type)]
ampli.info <- ampli.info[order(ampli.info$Type),]
for(cat in unique(ampli.info$Type)){
  sel.cols <- ampli.info$Type%in%cat
  sel.means <- apply(to.plot[,sel.cols,drop=FALSE],2,mean,na.rm=TRUE)
  to.plot[,sel.cols] <- to.plot[,sel.cols,drop=FALSE][,order(sel.means),drop=FALSE]
}
pheatmap(to.plot,
         color=rev(inferno(50)),
         annotation_col = ampli.info[,c('Type','GC_Content'),drop=FALSE],
         show_rownames=FALSE,show_colnames=FALSE,
         clustering_method = 'ward.D',
         cluster_cols = FALSE)
