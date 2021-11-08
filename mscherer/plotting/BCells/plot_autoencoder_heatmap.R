############## plot_autoencoder_heatmap.R ##############
#' With this script, we use the output of the methylation autoencoder, i.e., the mixture parameter of the forground and background read counts
#' distributions and plot it as a heatmap. The original read counts also have to be specified, as well as information about the amplicons.
#' An assignment of the cells to cell types also has to be specified and is used for cluster annotation.
#' Optionally, the matrix can be binarized to have an even cleaner picture.

library(pheatmap)
library(viridis)
sample <- 'BCells_Sample7_70_percent_good_performance'
meth.data <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/methylation_autoencoder/mixture_prob_dca.csv'),
                       header = TRUE,
                       row.names = 1)
input <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                   sep='\t',
                   header=TRUE)
ampli.info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.tsv')
mis.genes <- read.csv((paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/methylation_autoencoder/missing_amplicons.csv')))
mis.genes <- as.character(mis.genes[,2])
input <- input[,!(colnames(input)%in%mis.genes)]
#ampli.info <- ampli.info[!(row.names(ampli.info)%in%mis.genes),]
ampli.info <- ampli.info[colnames(input),]
colnames(meth.data) <- colnames(input)
row.names(meth.data) <- row.names(input)
ampli.info$Type <- as.factor(ampli.info$Type)
meth.data <- meth.data[,row.names(ampli.info)]
meth.data <- meth.data[,order(ampli.info$Type)]
ampli.info <- ampli.info[order(ampli.info$Type),]
for(cat in unique(ampli.info$Type)){
  sel.cols <- ampli.info$Type%in%cat
  sel.means <- apply(meth.data[,sel.cols,drop=FALSE],2,mean,na.rm=TRUE)
  meth.data[,sel.cols] <- meth.data[,sel.cols,drop=FALSE][,order(sel.means),drop=FALSE]
}
pheatmap(meth.data,
         color=rev(inferno(50)),
         annotation_col = ampli.info[,c('Type','GC_Content'),drop=FALSE],
         show_rownames=FALSE,show_colnames=FALSE,
         clustering_method = 'ward.D',
         cluster_cols = FALSE)
