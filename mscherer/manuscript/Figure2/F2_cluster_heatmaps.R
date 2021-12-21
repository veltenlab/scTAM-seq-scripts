################### F2_cluster_heatmaps.R ################################
#' This file generates heatmaps of clusters and associated differential CpGs.

library(pheatmap)
library(viridis)
sample <- 'Sample7_70_percent_good_performance_HP'
color_map <- list(CellType=c('naive'='#fcbd7e',
                             'memory1'='#fc6571',
                             'memory2'='#fc3262'),
                  Bulk=c('Naive B-cell high'='#fcbdbd',
                         'Memory B-cell high'='#fc00bd'),
                  Class=c('Naive B-cell high'='#fcbdbd',
                          'Memory B-cell high'='#fc00bd'),
                  CellType_reclustering=c('Cluster1'='#fcbd7e',
                                          'Cluster2a'='#fc3262',
                                          'Cluster2b'='#bd0071',
                                          'Cluster2c'='#8e008e'))
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/cluster_heatmaps/'
rowinfo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/rowinfo.csv'),
                    row.names = 1)
filtered.counts <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/",
  sample, '/tsv/', sample, '.barcode.cell.distribution.tsv'),
  row.names = 1,
  header=T)
filtered.counts <- ifelse(filtered.counts>0, 1, 0)
load('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/meth.data.numeric.Rdata')
cell_assignment <- read.table('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/MBC_assignment.txt')
meth.data.numeric <- meth.data.numeric[,as.character(cell_assignment$V5)]
for(clust1 in unique(rowinfo$CellType_reclustering)){
  for(clust2 in unique(rowinfo$CellType_reclustering)){
    diff_cpg_file <- paste0('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/differential/differential_CpGs_', clust1, 'vs', clust2, '.csv')
    if(!file.exists(diff_cpg_file)) next
    diff_cpgs <- read.csv(diff_cpg_file)
    mean_classes <- aggregate(t(meth.data.numeric[diff_cpgs$CpGID, ]), by=list(cell_assignment$V2), mean)
    row.names(mean_classes) <- mean_classes$Group.1
    mean_classes <- mean_classes[,-1]
    colnames(mean_classes) <- diff_cpgs$Amplicon
    png(paste0(plot_path, 'heatmap_', clust1, 'vs', clust2, '.png'),
        width=1000,
        height=1600)
    pheatmap(filtered.counts[row.names(rowinfo[rowinfo$CellType_reclustering%in%c(clust1, clust2), ]), diff_cpgs$Amplicon],
             annotation_row = subset(rowinfo, select=c('CellType_reclustering', 'Nfeatures')),
             annotation_col = as.data.frame(t(mean_classes)), 
             show_rownames = FALSE,
             show_colnames = FALSE,
             color=rev(inferno(50)),
             annotation_colors = color_map,
             fontsize=15)
    dev.off()
  }
}
