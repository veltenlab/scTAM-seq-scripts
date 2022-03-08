################### identify_doublets_per_cluster.R ################################
#' This file is used to identify further doublets that show substantial deviation from
#' the expected number of features based on a previous clustering of cells.

library(ComplexHeatmap)
library(viridis)
library(ggplot2)
#library(autothresholdr)
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
sample <- 'Sample11_70_percent_good_performance'
uncut <- 'Sample12_70_percent_good_performance'
protein <- 'Sample11'
plot_path <- '/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure2/Sample11/'

selected_amplicons <- read.table(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/',
                                       uncut, '/tsv/selected_amplicons.tsv'),
                              row.names = 1)
#selected_amplicons <- read.table('/users/lvelten/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv',
#                                  row.names = 1)
#selected_amplicons <- subset(selected_amplicons, subset = Type.of.amplicon=="CpG.B.cell.diff")
#rowinfo <- read.table(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/', sample, '/tsv/all_doublets.tsv'),
#                    row.names = 1)
#rowinfo <- read.csv(paste0('/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/all_doublets.csv'),
#                    row.names = 1)
rowinfo <- read.csv(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/', sample, '/tsv/rowinfo.csv'),
                     row.names = 1)
#rowinfo <- read.csv(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/', sample, '/tsv/doublet_scores_DoubletDetection.csv'),
#                    row.names = 2)
rowinfo <- rowinfo[rowinfo$Doublet==0, ]
filtered.counts <- read.table(paste0("/users/lvelten/project/Methylome/analysis/missionbio/BM/", sample, "/tsv/", sample, ".barcode.cell.distribution.tsv"),
                                     row.names = 1, header=T)

bulk.methylation <- read.table("/users/lvelten/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt", header=T)
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
#selected_data <- ifelse(filtered.counts[, row.names(selected_amplicons)]>0, 1, 0)
#selected_data <- ifelse(filtered.counts[row.names(rowinfo), ]>0, 1, 0)
#selected_data <- selected_data[which(rowinfo$Doublet==0), ]
#rowinfo <- rowinfo[which(rowinfo$Doublet==0), ]
# re_clustering <- cutree(hclust(dist(selected_data, 'binary'), method='ward.D2'), 4)
# clust_ass <- c('1'='Cluster4',
#                '2'='Cluster1',
#                '3'='Cluster2',
#                '4'='Cluster3')
# rowinfo[names(re_clustering), 'Cluster'] <- clust_ass[re_clustering]
# clust_ass <- c('1'='memory B-cells',
#                '2'='S2-S4 cells',
#                '3'='naive B-cells',
#                '4'='S1 cells')
# rowinfo[names(re_clustering), 'CellType'] <- clust_ass[re_clustering]
#' clust_ass <- c('1'='ns-memory B-cells',
#'                '2'='naive B-cells',
#'                '3'='cs-memory B-cells')#,
#'                #'4'='naive B-cells')
#' rowinfo[names(re_clustering), 'CellType'] <- clust_ass[re_clustering]
#rowinfo$Nfeatures <- apply(selected_data, 1, function(x){
#  sum(x>0)
#})
# prot_data <- read.table(paste0('/users/lvelten/project/Methylome/analysis/missionbio/tapestri/protein/',
#    protein, '/pipeline/results/tsv/Sample11-counts.tsv'),
#    header=TRUE)
# prot_data <- prot_data[prot_data$cell_barcode%in%row.names(rowinfo), ]
# prot_data <- split(prot_data, prot_data$ab_description)
# for(i in 1:length(prot_data)){
#   x <- prot_data[[i]]
#   row.names(x) <- x$cell_barcode
#   ab_counts <- x[row.names(rowinfo), 'raw']
#   rowinfo[, unique(x[, 'ab_description'])] <- ab_counts
# }
anno_col <- HeatmapAnnotation(subset(colinfo, select = c("S1.mean",
                                                         "S2.mean",
                                                         #"S3.mean",
                                                         "S4.mean",
                                                         #"GC.mean",
                                                         #"PB.mean",
                                                         "NBC.mean",
                                                         "MBC.mean")))
anno_row <- HeatmapAnnotation(subset(rowinfo, select = c("CellType_detailed")),
                              #"CD10",
                              #"CD27")),
                              #"CD19", 
                              #"CD34",
                              #"CD45RA")))
png(file.path(plot_path, 'heatmap_complete.png'),
    width=1000,
    height=1600)
Heatmap(selected_data, 
               annotation_col = anno_col,
               annotation_row = anno_row, 
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
ph <- pheatmap(selected_data, 
               annotation_col = subset(colinfo, select = c("S1.mean",
                                                           "S2.mean",
                                                           "S3.mean",
                                                           "S4.mean",
                                                           #"GC.mean",
                                                           #"PB.mean",
                                                           "NBC.mean",
                                                           "MBC.mean")),
               annotation_row = subset(rowinfo, select = c("CellType_detailed")),
               #"CD10",
               #"CD27")),
               #"CD19", 
               #"CD34",
               #"CD45RA")), 
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
s2_s4_clust <- hclust(dist(selected_data, 'binary'), 'ward.D2')
ct_detailed <- rowinfo$CellType
names(ct_detailed) <- row.names(rowinfo)
s2_s4_clust <- cutree(s2_s4_clust, 2)
n_clust <- names(s2_s4_clust)
s2_s4_clust<- c('1'='S2 cells',
  '2'='S3-S4 cells')[s2_s4_clust]
names(s2_s4_clust) <- n_clust
ct_detailed[names(s2_s4_clust)] <- s2_s4_clust
rowinfo[names(ct_detailed), 'CellType_detailed'] <- ct_detailed
#write.csv(rowinfo, paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/', sample, '/tsv/rowinfo.csv'))
#write.csv(colinfo, paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/', sample, '/tsv/colinfo.csv'))

to_plot <- rowinfo[, -c(1:6, 8)]
to_plot <- reshape2::melt(to_plot, id='CellType')
plot <- ggplot(to_plot, aes(x=CellType, y=value))+geom_boxplot()+facet_wrap(variable~., nrow=5)+plot_theme
ggsave(paste0(plot_path, 'AB_boxplots_cell_types.pdf'), plot)
apply(rowinfo[, -c(1:8)], 2, function(x){
  kruskal.test(x, g = rowinfo$CellType)$p.value
})
to_plot <- rowinfo[, -c(1:4, 6:8)]
to_plot <- reshape2::melt(to_plot, id='Cluster')
plot <- ggplot(to_plot, aes(x=Cluster, y=value))+geom_boxplot()+facet_wrap(variable~., nrow=5)+plot_theme
ggsave(paste0(plot_path, 'AB_boxplots_cluster.pdf'), plot)
apply(rowinfo[, -c(1:8)], 2, function(x){
  kruskal.test(x, g = rowinfo$Cluster)$p.value
})

res <- apply(selected_data, 1, function(x){
   sum(x>0)
})
