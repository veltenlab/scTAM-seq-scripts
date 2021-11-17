################### F1_F_all_scatterplots.R ################### 
#' This file generates pairwise scatterplots of bulk and single-cell data.

library(ggplot2)
library(GGally)
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=15),
                    axis.text=element_text(color='black',size=15),
                    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=12),
                    axis.text.y=element_text(size=12),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    legend.position='none')
color_map <- c('Cluster1'='#fcbd7e',
                                     'Cluster2a'='#e37e71',
                                     'Cluster2b'='#7f47bd',
                                     'Cluster2c'='#bd3262')

filtered.counts <- read.table("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/tsv/BCells_Sample7_70_percent_good_performance.barcode.cell.distribution_with_MCL.tsv", row.names = 1, header=T)
cell_metadata <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/rowinfo_reclustering_doublet.csv',
                          row.names=1)
amplicon.info <- read.csv("/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/selected_amplicons.csv", row.names = 1)

selected <- ifelse(filtered.counts[row.names(cell_metadata), row.names(amplicon.info)]>0, 1, 0)

load('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/meth.data.numeric.Rdata')
cell_assignment <- read.table('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/MBC_assignment.txt')
meth.data.numeric <- data.frame(meth.data.numeric[, grepl('NBC.PB', colnames(meth.data.numeric))],
                                meth.data.numeric[,as.character(cell_assignment$V5)])
more_info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
row.names(more_info) <- more_info$amplicon
joint_names <- intersect(row.names(more_info), colnames(selected))
selected <- selected[,joint_names]
colnames(selected) <- more_info[colnames(selected), 'background.cpgs']
selected <- selected[, !is.na(colnames(selected))]
mean_classes <- aggregate(t(meth.data.numeric[colnames(selected), ]), by=list(c(rep('naiveB', 5),cell_assignment$V2)), mean)
mean_cluster1 <- apply(selected[row.names(rowinfo)[rowinfo$CellType_reclustering%in%'Cluster1'], ], 2, function(x){
  sum(x>0)/length(x)
})
mean_cluster2a <- apply(selected[row.names(rowinfo)[rowinfo$CellType_reclustering%in%'Cluster2a'], ], 2, function(x){
  sum(x>0)/length(x)
})
mean_cluster2b <- apply(selected[row.names(rowinfo)[rowinfo$CellType_reclustering%in%'Cluster2b'], ], 2, function(x){
  sum(x>0)/length(x)
})
mean_cluster2c <- apply(selected[row.names(rowinfo)[rowinfo$CellType_reclustering%in%'Cluster2c'], ], 2, function(x){
  sum(x>0)/length(x)
})
mean_classes <- t(mean_classes)
colnames(mean_classes) <- mean_classes[1, ]
mean_classes <- as.data.frame(mean_classes[-1, ])
mean_classes$csMBC <- as.numeric(mean_classes$csMBC)
mean_classes$ncsMBC <- as.numeric(mean_classes$ncsMBC)
mean_classes$naiveB <- as.numeric(mean_classes$naiveB)
to_plot <- data.frame(mean_classes, 
                      'Cluster1'=mean_cluster1,
                      'Cluster2a'=mean_cluster2a,
                      'Cluster2b'=mean_cluster2b,
                      'Cluster2c'=mean_cluster2c)
plot <- ggpairs(to_plot)+plot_theme
ggsave(file.path(plot.path, 'all_scatterplots.png'), plot, width=200, height=150, units='mm')