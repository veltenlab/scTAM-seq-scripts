################### F1_A_heatmap.R ################################
#' This file is used to identify further doublets that show substantial deviation from
#' the expected number of features based on a previous clustering of cells.

library(pheatmap)
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
sample <- 'Sample11_70_percent_good_performance'
uncut <- 'Sample12_70_percent_good_performance'
protein <- 'Sample11'
plot_path <- '~'

selected_amplicons <- read.table(paste0('../../misc/', uncut, '/tsv/selected_amplicons.tsv'),
                              row.names = 1)
rowinfo <- read.csv(paste0('../../misc', sample, '/tsv/rowinfo.csv'),
                     row.names = 1)
rowinfo <- rowinfo[rowinfo$Doublet==0, ]
filtered.counts <- read.table(paste0("../../data/", sample, "/tsv/", sample, ".barcode.cell.distribution.tsv"),
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
ph <- pheatmap(selected_data, 
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
               color=rev(inferno(50)),
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
                                                           "NBC.mean",
                                                           "MBC.mean")),
               annotation_row = subset(rowinfo, select = c("Nfeatures",
                                                           "CellType_detailed")),
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
