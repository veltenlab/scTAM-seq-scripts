################### identify_doublets_per_cluster.R ################################
#' This file is used to identify further doublets that show substantial deviation from
#' the expected number of features based on a previous clustering of cells.

library(pheatmap)
library(viridis)
library(rtracklayer)
color_map <- list(CellType=c('naive'='#fcbd7e',
                             'memory1'='#fc6571',
                             'memory2'='#fc3262'),
                  Bulk=c('Naive B-cell high'='#fcbdbd',
                         'Memory B-cell high'='#fc00bd'),
                  Class=c('Naive B-cell high'='#fcbdbd',
                          'Memory B-cell high'='#fc00bd'),
                  CellType_reclustering=c('Cluster1'='#fcbd7e',
                             'Cluster2a'='#fc3262',
                             'Cluster2b'='#fc3262',
                             'Cluster2c'='#fc3262'))
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/'

rowinfo <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/rowinfo_reclustering_doublet.csv',
                    row.names = 1)
cell.info <- read.csv('/users/mscherer/cluster/project/Methylome/infos/BCells/cell_sexes_Sample7.csv',row.names = 1)
rowinfo$Sex <- cell.info[row.names(rowinfo), 'sex']
selected_amplicons <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/selected_amplicons.csv',
                               row.names = 1)
filtered.counts <- read.table("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/tsv/BCells_Sample7_70_percent_good_performance.barcode.cell.distribution_with_MCL.tsv", row.names = 1, header=T)
#amplicon.info <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv", header=T, row.names = 1)
#selected_amplicons <- amplicon.info[amplicon.info$Type.of.amplicon%in%"CpG.B.cell.diff", ]

bulk.methylation <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt", header=T)
row.names(bulk.methylation) <- bulk.methylation$amplicon
colinfo <- bulk.methylation[row.names(selected_amplicons), c('NBC.mean', 'MBC.mean', 'S1.mean', 'S2.mean', 'S3.mean', 'S4.mean')]
colinfo$Class <- ifelse(colinfo$NBC.mean>colinfo$MBC.mean, 'Naive B-cell high', 'Memory B-cell high')
atac_naive_clust_1 <- import('/users/mscherer/cluster/project/Methylome/data/external/ATAC/Satpathy/Naive_B_I-TileSize-500-normMethod-ReadsInTSS-ArchR.bw')
atac_naive_clust_2 <- import('/users/mscherer/cluster/project/Methylome/data/external/ATAC/Satpathy/Naive_B_II-TileSize-500-normMethod-ReadsInTSS-ArchR.bw')
atac_naive_clust_3 <- import('/users/mscherer/cluster/project/Methylome/data/external/ATAC/Satpathy/Naive_B_III-TileSize-500-normMethod-ReadsInTSS-ArchR.bw')
more_info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt',
                        row.names=8)
panel_gr <- makeGRangesFromDataFrame(more_info[row.names(colinfo), ],
                                      seqnames.field='Chr.hg19',
                                      start.field='Start.hg19',
                                      end.field='End.hg19')
op_1 <- findOverlaps(panel_gr, atac_naive_clust_1)
colinfo$ATAC_naiveB_I <- log(values(atac_naive_clust_1)[subjectHits(op_1), 'score']+1)
op_2 <- findOverlaps(panel_gr, atac_naive_clust_2)
colinfo$ATAC_naiveB_II <- log(values(atac_naive_clust_2)[subjectHits(op_2), 'score']+1)
op_3 <- findOverlaps(panel_gr, atac_naive_clust_3)
colinfo$ATAC_naiveB_III <- log(values(atac_naive_clust_3)[subjectHits(op_3), 'score']+1)
diff_cpgs <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/differential_naive/differential_CpGs_Cluster1avsCluster1b.csv')
colinfo$Differential <- ifelse(row.names(colinfo)%in%diff_cpgs$X, 'differential', 'other')
selected_data <- ifelse(filtered.counts[row.names(rowinfo), row.names(selected_amplicons)]>0, 1, 0)
#selected_data <- ifelse(filtered.counts[row.names(rowinfo), ]>0, 1, 0)
selected_data <- selected_data[rowinfo$Doublet%in%'Singlet', ]
re_clustering <- cutree(hclust(dist(selected_data, 'binary'), method='ward.D2'), 4)
clust_ass <- c('1'='Cluster2c',
               '2'='Cluster1',
               '3'='Cluster2b',
               '4'='Cluster2a')
#rowinfo[names(re_clustering), 'CellType_reclustering'] <- clust_ass[re_clustering]
png(file.path(plot_path, 'heatmap_complete_reclustering_naive_only.png'),
    width=1000,
    height=1600)
ph <- pheatmap(selected_data[rowinfo$CellType_reclustering%in%'Cluster1', ], 
               annotation_col = subset(colinfo, select = c('NBC.mean', 'MBC.mean', 'S1.mean', 'S2.mean', 'S3.mean', 'S4.mean', 'ATAC_naiveB_I', 'ATAC_naiveB_II', 'ATAC_naiveB_III', 'Differential')),
               annotation_row = subset(rowinfo, select = c("CellType_reclustering", "Nfeatures")), 
               clustering_distance_cols = "binary", 
               clustering_distance_rows = "binary", 
               show_rownames = F, 
               show_colnames = F, 
               cutree_rows = 2, 
               clustering_method = "complete",
               color=rev(inferno(50)),
               annotation_colors = color_map,
               fontsize=15)
dev.off()
#write.csv(rowinfo, '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/rowinfo_reclustering_doublet.csv')
clust <- hclust(dist(selected_data[rowinfo$CellType_reclustering%in%'Cluster1', ], 'binary'))
clust <- c('1'='Cluster1a',
           '2'='Cluster1b')[cutree(clust, 2)]
rowinfo <- rowinfo[rowinfo$CellType_reclustering%in%'Cluster1', ]
rowinfo$NaiveClustering <- clust
write.csv(rowinfo, '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/rowinfo_naive_clustering.csv')
