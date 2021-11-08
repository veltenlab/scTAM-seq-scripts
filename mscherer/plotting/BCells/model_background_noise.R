.libPaths(c(.libPaths(), '/users/mscherer/conda/envs/rnbeads/lib/R/library/'))
require(ggplot2)
require(plyr)
require(reshape2)
library(pheatmap)
filtered.counts <- read.table("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/tsv/BCells_Sample7_70_percent_good_performance.barcode.cell.distribution_with_MCL.tsv", row.names = 1, header=T)
amplicon.info <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv", header=T, row.names = 1)


bulk.methylation <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt", header=T)

autoencoder.bottleneck <- read.csv("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/methylation_autoencoder/bottleneck_good.csv", row.names = 1)
autoencoder.output <- read.csv("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/methylation_autoencoder/mixture_prob_dca_with_MCL.csv", row.names = 1)
missing.amplicons <- read.csv("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/methylation_autoencoder/missing_amplicons_with_MCL.csv", row.names = 1)
filtered.counts <- filtered.counts[,!colnames(filtered.counts) %in% missing.amplicons$X0]
dimnames(autoencoder.output) <- dimnames(filtered.counts)

control_amplicons <- filtered.counts[, row.names(amplicon.info[amplicon.info$Type.of.amplicon%in%'NonHhaI',])]
control_amplicons_binary <- ifelse(control_amplicons>0, 1, 0)
pheatmap(log10(control_amplicons+1))
pheatmap(control_amplicons_binary, clustering_method = 'ward.D')
clusters <- cutree(hclust(dist(t(control_amplicons_binary)), method='ward.D'),2)
heterogeneous_cluster_amplicons <- names(which(clusters==1))
meth_cells <- apply(control_amplicons_binary[, heterogeneous_cluster_amplicons], 2, function(x){
  sum(x>0)/length(x)
})
hist(meth_cells)

memory_selected <- filtered.counts[row.names(rowinfo[rowinfo$Cluster%in%c(1),]),]
memory_selected <- memory_selected[,  joint_names]
memory_selected <- apply(memory_selected, 2, function(x){
  sum(x>0)/length(x)
})
hist(memory_selected)
