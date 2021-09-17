library(pheatmap)
library(viridis)
sample <- 'Sample5_80_percent'
meth.data <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/error_mod/pyro/first_try.csv',
                      header = TRUE,
                      row.names = 1)
input <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/Sample5_80_percent/tsv/Sample5_80_percent.barcode.cell.distribution.tsv'),
                    sep='\t',
                    header=TRUE)
ampli.info <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/amplicons/cell_lines/R/amplicon_info_with_RnBeads.csv')
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
row.names(meth.data) <- row.names(input)
clust.file <- clust.file[row.names(meth.data),]
colnames(meth.data) <- ampli.info$AmpID_design_1898
row.names(ampli.info) <- ampli.info$AmpID_design_1898
ampli.info$Type[grepl('Aci',ampli.info$Type)] <- 'Aci'
ampli.info$Type <- as.factor(ampli.info$Type)
meth.data <- meth.data[,row.names(ampli.info)]
meth.data <- meth.data[,order(ampli.info$Type)]
ampli.info <- ampli.info[order(ampli.info$Type),]
for(cat in unique(ampli.info$Type)){
  sel.cols <- ampli.info$Type%in%cat
  sel.means <- apply(meth.data[,sel.cols],2,mean,na.rm=TRUE)
  meth.data[,sel.cols] <- meth.data[,sel.cols][,order(sel.means)]
}
pheatmap(meth.data,
         annotation_row = clust.file[,'CellType',drop=FALSE],
         color=rev(inferno(50)),
         annotation_col = ampli.info[,c('Type','CG'),drop=FALSE],
         show_rownames=FALSE,show_colnames=FALSE,
         clustering_method = 'ward.D',
         cluster_cols = FALSE,
         annotation_colors=list(CellType=c('Jurkat'='#ff9289', 'K562'='#82b7ff', 'Mixed'='grey'),
                                Type=c('Aci'='grey','Hha Jurkat high'='#b3645eff','Hha K562 high'='#597daeff')))
