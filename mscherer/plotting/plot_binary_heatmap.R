############## plot_binary_heatmap.R ##############
#' This script is used to plot a binary heatmap of the read counts. A read cutoff has to be specified (default 10)
#' to differentiate between presence and absence of read counts, i.e., DNA methylation. Additionally, the read count
#' file and the assignment of cells to clusters has to be specified.

library(pheatmap)
cut.off <- 5
sample <- 'Sample5_80_percent'
dat.s3 <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'))
to.plot <- apply(dat.s3,2,function(x){
  ifelse(x>cut.off,1,0)
})
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
ampli.info <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/amplicons/cell_lines/R/amplicon_info_with_RnBeads.csv')
row.names(ampli.info) <- ampli.info$AmpID_design_1898
colnames(to.plot) <- colnames(dat.s3)
row.names(to.plot) <- row.names(dat.s3)
clust.file <- clust.file[row.names(to.plot),]
row.names(clust.file) <- row.names(to.plot)
ampli.info$Type[grepl('Aci',ampli.info$Type)] <- 'Aci'
ampli.info$Type <- as.factor(ampli.info$Type)
to.plot <- to.plot[,row.names(ampli.info)]
to.plot <- to.plot[,order(ampli.info$Type)]
ampli.info <- ampli.info[order(ampli.info$Type),]
to.plot <- to.plot[!(clust.file$CellType%in%'Mixed'),]
for(cat in unique(ampli.info$Type)){
  sel.cols <- ampli.info$Type%in%cat
  sel.means <- apply(to.plot[,sel.cols],2,mean,na.rm=TRUE)
  to.plot[,sel.cols] <- to.plot[,sel.cols][,order(sel.means)]
}
pheatmap(to.plot,
         annotation_row = clust.file[,'CellType',drop=FALSE],
         color=rev(inferno(50)),
         annotation_col = ampli.info[,c('Type','CG'),drop=FALSE],
         show_rownames=FALSE,show_colnames=FALSE,
         clustering_method = 'ward.D',
         cluster_cols = FALSE,
         annotation_colors=list(CellType=c('Jurkat'='#ff9289', 'K562'='#82b7ff', 'Mixed'='grey'),
                                Type=c('Aci'='grey','Hha Jurkat high'='#b3645eff','Hha K562 high'='#597daeff')))

library(pheatmap)
cut.off <- 10
sample <- 'Sample4_80_percent'
dat.s3 <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'))
sel.amplis <- c('AMPL130616',
                'AMPL130659',
                'AMPL130742',
                'AMPL130770',
                'AMPL134231',
                'AMPL130961',
                'AMPL134313',
                'AMPL134314')
dat.s3 <- dat.s3[,sel.amplis]
to.plot <- apply(dat.s3,2,function(x){
  ifelse(x>cut.off,1,0)
})
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
ampli.info <- read.csv('/users/mscherer/cluster/project/Methylome/infos/cell_lines/analysis_Agostina_Apr21/Summary_methylation_values_amplicons_K562_Jurkat_v2.csv')
row.names(ampli.info) <- ampli.info$AmpID_design_1898
pheatmap(to.plot,
         cluster_rows = TRUE,
         annotation_row = clust.file[,'CellType',drop=FALSE],
         annotation_col = ampli.info[,'Label',drop=FALSE],
         show_colnames = TRUE,
         show_rownames = FALSE)
