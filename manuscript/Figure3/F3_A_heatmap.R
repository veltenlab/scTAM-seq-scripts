################### F1_A_heatmap.R ################################
#' This file is used to identify further doublets that show substantial deviation from
#' the expected number of features based on a previous clustering of cells.

library(ComplexHeatmap)
library(viridis)
library(ggplot2)
color_map <- list(CellType_broad=c('naive B-cells'='#fcbd7e',
                             'memory B-cells'='#fc6571',
                             'memory B-cells'='#fc3262'),
                  Bulk=c('Naive B-cell high'='#fcbdbd',
                         'Memory B-cell high'='#fc00bd'),
                  Class=c('Naive B-cell high'='#fcbdbd',
                          'Memory B-cell high'='#fc00bd'),
                  CellType=c('naive B-cells'='#fcbd7e',
                             'memory B-cells'='#fc6571',
                             'S2-S4 cells'='#fcef7e',
                             'S1 cells'='#bdfc7e'),
                  CellType_detailed=c('naive B-cells'='#fcbd7e',
                             'memory B-cells'='#fc6571',
                             'S2 cells'='#d7ef7e',
                             'S3-S4 cells'='#fcef7e',
                             'S1 cells'='#bdfc7e'),
                  MBC.mean=rev(c('#43a2ca',
                               '#7bccc4',
                               '#bae4bc',
                               '#f0f9e8')),
                  NBC.mean=rev(c('#8856a7',
                               '#8c96c6',
                               '#b3cde3',
                               '#edf8fb')),
                  S3.mean=rev(c('#e34a33',
                            '#fc8d59',
                            '#fdcc8a',
                            '#fef0d9')),
                  S4.mean=rev(c('#e6550d',
                            '#fd8d3c',
                            '#fdbe85',
                            '#feedde')),
                  S1.mean=rev(c('#636363',
                            '#969696',
                            '#cccccc',
                            '#f7f7f7')),
                  S2.mean=rev(c('#31a354',
                            '#74c476',
                            '#bae4b3',
                            '#edf8e9')))
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(angle=90, hjust=1, color='black',size=5),
                    axis.ticks=element_line(color='black', size=.25),
                    strip.background = element_blank(),
                    #strip.text.x = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')
sample <- 'GSM5935921_BM_HhaI'
uncut <- 'GSM5935923_BM_undigested'
plot_path <- '~'

selected_amplicons <- read.table(paste0('../../misc/', uncut, '/tsv/selected_amplicons.tsv'),
                              row.names = 1)
rowinfo <- read.csv(paste0('../../misc/', sample, '/tsv/rowinfo.csv'),
                     row.names = 1)
rowinfo <- rowinfo[rowinfo$Doublet==0, ]
filtered.counts <- read.table(paste0("../../data/", sample, ".tsv.gz"),
                                     row.names = 1, header=T)

bulk.methylation <- read.table("../../misc/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt",
                               header=T)
row.names(bulk.methylation) <- bulk.methylation$amplicon
colinfo <- bulk.methylation[row.names(selected_amplicons), c("S1.mean",
                                                             "S2.mean",
                                                             "S3.mean",
                                                             "S4.mean",
                                                             "GC.mean",
                                                             "PB.mean",
                                                             "NBC.mean",
                                                             "MBC.mean")]
selected_data <- ifelse(filtered.counts[row.names(rowinfo), row.names(selected_amplicons)]>0, 1, 0)
png(file.path(plot_path, 'heatmap_complete.png'),
    width=1000,
    height=1600)
pheatmap(selected_data, 
               annotation_col = subset(colinfo, select = c("S1.mean",
                                                           "S2.mean",
                                                           "S3.mean",
                                                           "S4.mean",
                                                           "NBC.mean",
                                                           "MBC.mean")),
               annotation_row = subset(rowinfo, select = c("Nfeatures",
                                                           "CellType_detailed")),
               clustering_distance_cols = "binary", 
               clustering_distance_rows = "binary",
               show_rownames = F, 
               show_colnames = F, 
               cutree_rows = 4, 
               clustering_method = "ward.D2",
               col=rev(inferno(50)),
               annotation_colors = color_map,
               fontsize=15,
               legend = FALSE,
               annotation_legend = FALSE,
               labels_row = FALSE,
               labels_col = FALSE)
dev.off()

selected_data <- selected_data[which(rowinfo$CellType%in%'S2-S4 cells'), ]
png(file.path(plot_path, 'heatmap_S2_S4_cluster.png'),
    width=1000,
    height=1600)
pheatmap(selected_data, 
               annotation_col = subset(colinfo, select = c("S1.mean",
                                                           "S2.mean",
                                                           "S3.mean",
                                                           "S4.mean",
                                                           "NBC.mean",
                                                           "MBC.mean")),
               clustering_distance_cols = "binary", 
               clustering_distance_rows = "binary", 
               show_rownames = F, 
               show_colnames = F, 
               cutree_rows = 2, 
               clustering_method = "ward.D2",
               color=rev(inferno(50)),
               annotation_colors = color_map,
               fontsize=15,
               legend = FALSE,
               annotation_legend = FALSE,
               labels_row = FALSE,
               labels_col = FALSE)
dev.off()
