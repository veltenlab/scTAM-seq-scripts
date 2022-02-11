################### identify_doublets_per_cluster.R ################################
#' This file is used to identify further doublets that show substantial deviation from
#' the expected number of features based on a previous clustering of cells.

library(pheatmap)
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
                             'ns-memory B-cells'='#fc3262',
                             'cs-memory B-cells'='#8e008e'))
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(angle=90, hjust=1, color='black',size=5),
                    axis.ticks=element_line(color='black', size=.25),
                    strip.background = element_blank(),
                    #strip.text.x = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')
sample <- 'Sample8_70_percent_good_performance'
uncut <- 'Sample6_70_percent_good_performance'
protein <- 'Sample8'
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/Sample8/'

#selected_amplicons <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/',
#                                       uncut, '/tsv/selected_amplicons.tsv'),
#                              row.names = 1)
selected_amplicons <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv',
                                  row.names = 1)
selected_amplicons <- subset(selected_amplicons, subset = Type.of.amplicon=="CpG.B.cell.diff")
#rowinfo <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/all_doublets.tsv'),
#                    row.names = 1)
#rowinfo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/all_doublets.csv'),
#                    row.names = 1)
rowinfo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/rowinfo.csv'),
                    row.names = 1)
filtered.counts <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/", sample, "/tsv/", sample, ".barcode.cell.distribution.tsv"),
                                     row.names = 1, header=T)

bulk.methylation <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt", header=T)
row.names(bulk.methylation) <- bulk.methylation$amplicon
colinfo <- bulk.methylation[row.names(selected_amplicons), c('NBC.mean', 'MBC.mean')]
colinfo$Class <- ifelse(colinfo$NBC.mean>colinfo$MBC.mean, 'Naive B-cell high', 'Memory B-cell high')
selected_data <- ifelse(filtered.counts[row.names(rowinfo), row.names(selected_amplicons)]>0, 1, 0)
#selected_data <- ifelse(filtered.counts[, row.names(selected_amplicons)]>0, 1, 0)
#selected_data <- ifelse(filtered.counts[row.names(rowinfo), ]>0, 1, 0)
#selected_data <- selected_data[which(rowinfo$Doublet==0), ]
#rowinfo <- rowinfo[which(rowinfo$Doublet==0), ]
#' re_clustering <- cutree(hclust(dist(selected_data, 'binary'), method='ward.D2'), 3)
#' clust_ass <- c('1'='Cluster2a',
#'                '2'='Cluster1',
#'                '3'='Cluster2b')
#'                #'4'='Cluster1')
#' rowinfo[names(re_clustering), 'Cluster'] <- clust_ass[re_clustering]
#' clust_ass <- c('1'='memory B-cells',
#'                '2'='naive B-cells',
#'                '3'='memory B-cells')#,
#'                #'4'='naive B-cells')
#' rowinfo[names(re_clustering), 'CellType_broad'] <- clust_ass[re_clustering]
#' clust_ass <- c('1'='ns-memory B-cells',
#'                '2'='naive B-cells',
#'                '3'='cs-memory B-cells')#,
#'                #'4'='naive B-cells')
#' rowinfo[names(re_clustering), 'CellType'] <- clust_ass[re_clustering]
#' rowinfo$Nfeatures <- apply(selected_data, 1, function(x){
#'   sum(x>0)
#' })
# prot_data <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/protein/',
#   protein, '/pipeline/results/tsv/Sample8-counts.tsv'),
#   header=TRUE)
# prot_data <- prot_data[prot_data$cell_barcode%in%row.names(rowinfo), ]
# prot_data <- split(prot_data, prot_data$ab_description)
# for(i in 1:length(prot_data)){
#   x <- prot_data[[i]]
#   row.names(x) <- x$cell_barcode
#   ab_counts <- x[row.names(rowinfo), 'raw']
#   rowinfo[, unique(x[, 'ab_description'])] <- ab_counts
#   thres.hold <- auto_thresh(ab_counts, method='otsu', ignore_na = TRUE)
#   rowinfo[, paste0(unique(x[, 'ab_description']), '_binary')] <- ifelse(ab_counts>thres.hold,
#                                                                         paste0(unique(x[, 'ab_description']), '+'),
#                                                                         paste0(unique(x[, 'ab_description']), '-')
#   )
# }
png(file.path(plot_path, 'heatmap_complete.png'),
    width=1000,
    height=1600)
ph <- pheatmap(selected_data, 
               annotation_col = subset(colinfo, select = c("NBC.mean",
                                                           "MBC.mean")),
               annotation_row = subset(rowinfo, select = c("CellType",
                                                           "Nfeatures")),
                                                           #"CD10",
                                                           #"CD27")),
                                                           #"CD19", 
                                                           #"CD34",
                                                           #"CD45RA")), 
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
#write.csv(rowinfo, paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/rowinfo.csv'))
#write.csv(colinfo, paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/colinfo.csv'))

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
