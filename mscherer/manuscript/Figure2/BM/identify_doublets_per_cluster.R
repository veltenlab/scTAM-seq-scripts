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
plot_path <- '/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure2/Sample11/'
sample <- 'Sample11_70_percent_good_performance'
uncut <- 'Sample12_70_percent_good_performance'

rowinfo <- read.csv(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/', sample, '/tsv/rowinfo.csv'),
                  row.names = 1)
#rowinfo <- read.table(paste0('/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/rowinfo.tsv'),
#                    row.names = 1)
#rowinfo <- read.csv(paste0('/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv//doublet_scores_DoubletDetection.csv'),
#                                      row.names = 2)
#rowinfo <- rowinfo[which(rowinfo$DoubletDetectionLabel==0), ]
selected_amplicons <- read.table(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/', 
                                        uncut, '/tsv/selected_amplicons.tsv'),
                                 row.names = 1)
filtered.counts <- read.table(paste0("/users/lvelten/project/Methylome/analysis/missionbio/BM/", 
  sample, "/tsv/", sample, ".barcode.cell.distribution.tsv"),
  row.names = 1,
  header=T)

selected_data <- ifelse(filtered.counts[row.names(rowinfo), row.names(selected_amplicons)]>0, 1, 0)

ph <- pheatmap(selected_data, 
        #       annotation_col = subset(fpr, select = c("Bulk")),
               annotation_row = subset(rowinfo, select = c("DoubletDetectionLabel")), 
               clustering_distance_cols = "binary", 
               clustering_distance_rows = "binary", 
               show_rownames = F, 
               show_colnames = F, 
               cutree_rows = 4, 
               clustering_method = "ward.D2",
               color=rev(inferno(50)),
               annotation_colors = color_map,
               fontsize=15)

sel_dat <- subset(rowinfo, subset=Cluster=='Cluster1', select='Nfeatures')
hist(unlist(sel_dat))
pheatmap(selected_data[row.names(sel_dat), ], 
         #       annotation_col = subset(fpr, select = c("Bulk")),
         annotation_row = subset(rowinfo, select = c("DoubletDetectionLabel", "Nfeatures")), 
         cutree_rows = 4, 
         clustering_method = "ward.D2",
         color=rev(inferno(50)),
         show_rownames = F, 
         show_colnames = F,
         annotation_colors = color_map,
         fontsize=15)
clust <- hclust(dist(selected_data[row.names(sel_dat), ], 'binary'), method='ward.D2')
clust <- cutree(clust, 4)
pot_doublets <- names(clust[clust==4])

rowinfo$DoubletClustering <- rowinfo$DoubletDetectionLabel
rowinfo <- read.csv('/users/lvelten/project/Methylome/analysis/missionbio/BM/Sample11_70_percent_good_performance/tsv/doublet_scores_DoubletDetection.csv', row.names=2)
rowinfo$DoubletClustering <- rowinfo$DoubletDetectionLabel
rowinfo[pot_doublets, 'DoubletClustering'] <- 1
rowinfo$Doublet <- rowinfo$DoubletClustering
rowinfo <- rowinfo[, -1]
png('/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Supplement/heatmap_Sample11_doublet_marked.png',
    width=1000,
    height=1600)
pheatmap(selected_data, 
         #       annotation_col = subset(fpr, select = c("Bulk")),
         annotation_row = subset(rowinfo, select = c("DoubletDetectionLabel", "DoubletClustering")), 
         cutree_rows = 4, 
         clustering_method = "ward.D2",
         color=rev(inferno(50)),
         show_rownames = F, 
         show_colnames = F,
         annotation_colors = color_map,
         fontsize=15)
dev.off()

write.table(rowinfo,
          paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/',
                 sample, '/tsv/all_doublets.tsv'))
