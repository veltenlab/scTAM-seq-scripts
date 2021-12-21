################### identify_doublets_per_cluster.R ################################
#' This file is used to identify further doublets that show substantial deviation from
#' the expected number of features based on a previous clustering of cells.

library(pheatmap)
color_map <- list(CellType_broad=c('naive B-cells'='#fcbd7e',
                                   'memory B-cells'='#fc6571',
                                   'memory B-cells'='#fc3262'),
                  Bulk=c('Naive B-cell high'='#fcbdbd',
                         'Memory B-cell high'='#fc00bd'),
                  Class=c('Naive B-cell high'='#fcbdbd',
                          'Memory B-cell high'='#fc00bd'),
                  CellType=c('naive B-cells'='#fcbd7e',
                             'ns-memory B-cells'='#fc3262',
                             'ns-memory B-cells'='#fc3262',
                             'cs-memory B-cells'='#fc3262'))
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/Sample8/'
sample <- 'Sample8_70_percent_good_performance'
uncut <- 'Sample6_70_percent_good_performance'

rowinfo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/rowinfo.csv'),
                  row.names = 1)
selected_amplicons <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', 
                                        uncut, '/tsv/selected_amplicons.tsv'),
                                 row.names = 1)
filtered.counts <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/", 
  sample, "/tsv/", sample, ".barcode.cell.distribution.tsv"),
  row.names = 1,
  header=T)

selected_data <- ifelse(filtered.counts[row.names(rowinfo), row.names(selected_amplicons)]>0, 1, 0)

ph <- pheatmap(selected_data, 
        #       annotation_col = subset(fpr, select = c("Bulk")),
               annotation_row = subset(rowinfo, select = c("CellType", "Nfeatures")), 
               clustering_distance_cols = "binary", 
               clustering_distance_rows = "binary", 
               show_rownames = F, 
               show_colnames = F, 
               cutree_rows = 3, 
               clustering_method = "ward.D2",
               color=rev(inferno(50)),
               annotation_colors = color_map,
               fontsize=15)

sel_dat <- subset(rowinfo, subset=Cluster=='Cluster1', select='Nfeatures')
hist(unlist(sel_dat))
pot_doublets <- row.names(sel_dat)[sel_dat$Nfeatures>400]
sel_dat <- subset(rowinfo, subset=Cluster=='Cluster2a', select='Nfeatures')
hist(unlist(sel_dat))
pot_doublets <- row.names(sel_dat)[sel_dat$Nfeatures>180]
#pot_doublets <- c(row.names(sel_dat)[sel_dat$Nfeatures>180], pot_doublets)
doublets <- rep('Singlet', nrow(rowinfo))
doublets[row.names(rowinfo)%in%pot_doublets] <- 'Doublet'
rowinfo$Doublet <- doublets
png(file.path(plot.path, 'heatmap_complete_doublets_marked.png'),
    width=1000,
    height=1600)
ph <- pheatmap(selected_data, 
        #       annotation_col = subset(selected_amplicons, select = c("Bulk")),
               annotation_row = subset(rowinfo, select = c("CellType_reclustering", "Nfeatures", "Doublet")), 
               clustering_distance_cols = "binary", 
               clustering_distance_rows = "binary", 
               show_rownames = F, 
               show_colnames = F, 
               cutree_rows = 3, 
               clustering_method = "ward.D2",
               color=rev(inferno(50)),
               annotation_colors = color_map,
               fontsize=15)
dev.off()

selected_data <- selected_data[rowinfo$Doublet%in%'Singlet', ]
re_clustering <- cutree(hclust(dist(selected_data, 'binary'), method='ward.D2'), 3)
clust_ass <- c('1'='cs-memory B-cells',
               '2'='naive B-cells',
               '3'='ns-memory B-cells')
color_map <- list(CellType=c('naive'='#fcbd7e',
                             'memory1'='#fc6571',
                             'memory2'='#fc3262'),
                  Bulk=c('Naive B-cell high'='#fcbdbd',
                         'Memory B-cell high'='#fc00bd'),
                  CellType_reclustering=c('Cluster1'='#fcbd7e',
                                          'Cluster2a'='#fc3262',
                                          'Cluster2b'='#fc3262',
                                          'Cluster2c'='#fc3262'))
rowinfo[names(re_clustering), 'CellType_reclustering'] <- clust_ass[re_clustering]
png(file.path(plot.path, 'heatmap_complete_reclustering.png'),
    width=1000,
    height=1600)
ph <- pheatmap(selected_data, 
            #   annotation_col = subset(selected_amplicons, select = c("Bulk")),
               annotation_row = subset(rowinfo, select = c("CellType_reclustering", "Nfeatures")), 
               clustering_distance_cols = "binary", 
               clustering_distance_rows = "binary", 
               show_rownames = F, 
               show_colnames = F, 
               cutree_rows = 3, 
               clustering_method = "ward.D2",
               color=rev(inferno(50)),
               annotation_colors = color_map,
               fontsize=15)
dev.off()
rowinfo$DoubletClustering <- rowinfo$DoubletDetectionLabel
rowinfo <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/tsv/doublet_scores_DoubletDetection.csv', row.names=2)
rowinfo$DoubletClustering <- rowinfo$DoubletDetectionLabel
rowinfo[pot_doublets, 'DoubletClustering'] <- 1
rowinfo$Doublet <- rowinfo$DoubletClustering
rowinfo <- rowinfo[, -1]
write.table(rowinfo,
          paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/',
                 sample, '/tsv/all_doublets.tsv'))
write.csv(pot_doublets, file.path(plot.path,  paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/',
                                                     sample, 'reclustering_doublet_names.csv')))