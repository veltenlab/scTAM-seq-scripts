################### identify_doublets_per_cluster.R ################################
#' This file is used to identify further doublets that show substantial deviation from
#' the expected number of features based on a previous clustering of cells.

library(pheatmap)
library(viridis)
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
sample <- 'Sample8_70_percent_good_performance'
sample2 <- 'Sample7_70_percent_good_performance'
uncut <- 'Sample6_70_percent_good_performance'
protein <- 'Sample8'
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/'

selected_amplicons <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', 
                                        uncut, '/tsv/selected_amplicons.tsv'),
                               row.names = 1)
rowinfo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/doublet_scores_DoubletDetection.csv'),
                    row.names = 1)
rowinfo2 <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample2, '/tsv/doublet_scores_DoubletDetection.csv'),
                    row.names = 1)
rowinfo2$Barcode <- row.names(rowinfo2)
rowinfo <- as.data.frame(rbind(rowinfo, rowinfo2[, c('Barcode', 'DoubletDetectionScore', 'DoubletDetectionLabel')]))
filtered.counts <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/", sample, "/tsv/", sample, ".barcode.cell.distribution.tsv"),
                                     row.names = 1, header=T)
barcodes <- row.names(filtered.counts)
filtered.counts.second <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/", sample2, "/tsv/", sample2, ".barcode.cell.distribution.tsv"),
                              row.names = 1, header=T)
barcodes <- c(barcodes, row.names(filtered.counts.second))
filtered.counts <- as.data.frame(rbind(filtered.counts, filtered.counts.second))

bulk.methylation <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt", header=T)
row.names(bulk.methylation) <- bulk.methylation$amplicon
colinfo <- bulk.methylation[row.names(selected_amplicons), c('NBC.mean', 'MBC.mean')]
colinfo$Class <- ifelse(colinfo$NBC.mean>colinfo$MBC.mean, 'Naive B-cell high', 'Memory B-cell high')
#selected_data <- ifelse(filtered.counts[row.names(rowinfo), row.names(selected_amplicons)]>0, 1, 0)
selected_data <- ifelse(filtered.counts[, row.names(selected_amplicons)]>0, 1, 0)
#selected_data <- ifelse(filtered.counts[row.names(rowinfo), ]>0, 1, 0)
# selected_data <- selected_data[which(rowinfo$Doublet==0), ]
# rowinfo <- rowinfo[which(rowinfo$Doublet==0), ]
colliding.barcodes <- plyr::count(rowinfo$Barcode)
colliding.barcodes <- colliding.barcodes[colliding.barcodes$freq>1, 'x']
rowinfo <- rowinfo[-(which(rowinfo$Barcode%in%colliding.barcodes)), ]
selected_data <- selected_data[-(which(barcodes%in%colliding.barcodes)), ]
barcodes <- barcodes[-(which(barcodes%in%colliding.barcodes))]
row.names(rowinfo) <- rowinfo$Barcode
row.names(selected_data) <- barcodes
selected_data <- selected_data[which(rowinfo$DoubletDetectionLabel==0), ]
rowinfo <- rowinfo[which(rowinfo$DoubletDetectionLabel==0), ]
re_clustering <- cutree(hclust(dist(selected_data, 'binary'), method='ward.D2'), 4)
clust_ass <- c('1'='Cluster2a',
               '2'='Cluster2b',
               '3'='Cluster2c',
               '4'='Cluster1')
rowinfo[names(re_clustering), 'Cluster'] <- clust_ass[re_clustering]
clust_ass <- c('1'='memory B-cells',
               '2'='memory B-cells',
               '3'='memory B-cells',
               '4'='naive B-cells')
rowinfo[names(re_clustering), 'CellType_broad'] <- clust_ass[re_clustering]
clust_ass <- c('1'='cs-memory B-cells',
               '2'='ns-memory B-cells',
               '3'='ns-memory B-cells',
               '4'='naive B-cells')
rowinfo[names(re_clustering), 'CellType'] <- clust_ass[re_clustering]
prot_data <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/protein/', 
                               protein, '/pipeline/results/tsv/Sample8-counts.tsv'),
                        header=TRUE)
prot_data <- prot_data[prot_data$cell_barcode%in%row.names(rowinfo), ]
prot_data <- split(prot_data, prot_data$ab_description)
for(i in 1:length(prot_data)){
  x <- prot_data[[i]]
  row.names(x) <- x$cell_barcode
  ab_counts <- log10(x[row.names(rowinfo), 'raw'])
  rowinfo[, unique(x[, 'ab_description'])] <- ab_counts
}
rowinfo$Nfeatures <- apply(selected_data, 1, function(x){
  sum(x>0)
})
png(file.path(plot_path, 'heatmap_complete_Sample8andSample7.png'),
    width=1000,
    height=1600)
ph <- pheatmap(selected_data, 
               annotation_col = colinfo,
               annotation_row = subset(rowinfo, select = c("CellType",
                                                           "Nfeatures",
                                                           "CD10",
                                                           "CD11c",
                                                           "CD25",
                                                           "CD27",
                                                           "CD38",
                                                           "CD62L")), 
               clustering_distance_cols = "binary", 
               clustering_distance_rows = "binary", 
               show_rownames = F, 
               show_colnames = F, 
               cutree_rows = 4, 
               clustering_method = "ward.D2",
               color=rev(inferno(50)),
               annotation_colors = color_map,
               fontsize=15)
dev.off()
write.csv(rowinfo, paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/rowinfo.csv'))
write.csv(colinfo, paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/colinfo.csv'))