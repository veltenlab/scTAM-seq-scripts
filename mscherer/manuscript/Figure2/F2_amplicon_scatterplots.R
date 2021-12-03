################### F2_amplicon_scatterplots.R ################### 
#' This file generates pairwise scatterplots of bulk and single-cell data.

library(ggplot2)
library(pheatmap)
library(viridis)
library(GGally)
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/amplicon_correlation/'
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
#amplicon.info <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt", row.names = 1)
#amplicon.info <- subset(amplicon.info, subset = Type.of.amplicon=='CpG.B.cell.diff')
#row.names(amplicon.info) <- amplicon.info$amplicon
# more_info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
# row.names(more_info) <- more_info$amplicon
# all.cpgs <- more_info[row.names(amplicon.info), 'background.cpgs']
# load('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/meth.data.numeric.Rdata')
# cell_assignment <- read.table('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/MBC_assignment.txt')
# meth.data.numeric <- meth.data.numeric[,as.character(cell_assignment$V5)]
# means_csMBC <- apply(meth.data.numeric[all.cpgs, cell_assignment$V2=='csMBC'], 1, mean)
# means_ncsMBC <- apply(meth.data.numeric[all.cpgs, cell_assignment$V2=='ncsMBC'], 1, mean)
# selected_amplicons <- (means_csMBC>0.25&means_csMBC<0.75)&(means_ncsMBC>0.25&means_ncsMBC<0.75)
# #selected_amplicons[sample(1:length(selected_amplicons), 47)] <- FALSE
# bulk_table <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
# means_bulk <- bulk_table[, c('NBC.mean', 'MBC.mean')]
# colnames(means_bulk) <- c('naiveB', 'memoryB')
# row.names(means_bulk) <- bulk_table$amplicon
# means_bulk <- means_bulk[row.names(amplicon.info), ]
# selected_amplicons <- (means_bulk$naiveB>0.25&means_bulk$naiveB<0.75)&(means_bulk$memoryB>0.25&means_bulk$memoryB<0.75)
cluster_meth <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/dropou_modeling/all_corrected_stan_clusters.csv', row.names=1)
selected <- ifelse(filtered.counts[row.names(cell_metadata), row.names(amplicon.info)]>0, 1, 0)
cluster_meth <- cluster_meth[row.names(amplicon.info), ]
selected_data <- selected_data[, row.names(amplicon.info)]

for(clust in unique(cell_metadata$CellType_reclustering)){
  selected_amplicons <- cluster_meth[, clust]>0.25&cluster_meth[, clust]<0.75
#  vars <- apply(selected[cell_metadata$CellType_reclustering%in%clust, ], 2, var)
#  selected_amplicons <- rank(vars)>(length(vars)-30)|rank(vars)<10
#  selected_amplicons <- rank(vars)>(length(vars)-20)#|rank(vars)<10
  to_write <- data.frame(Amplicon=row.names(amplicon.info)[selected_amplicons], CpGID=more_info[row.names(amplicon.info)[selected_amplicons], 'background.cpgs'])
  write.csv(to_write,
            paste0('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/amplicon_correlation/variable_amplicons_', clust, '.csv'))
  amplicon.info$Selected <- ifelse(selected_amplicons, 'Yes', 'No')
  
  png(file.path(plot_path, paste0('heatmap_selected_marked_', clust, '.png')))
  pheatmap(selected_data, 
                 annotation_col = amplicon.info[, c('Selected'), drop=FALSE],
                 clustering_distance_cols = "binary", 
                 clustering_distance_rows = "binary", 
                 show_rownames = F, 
                 show_colnames = F, 
                 cutree_rows = 4, 
                 clustering_method = "ward.D2",
                 color=rev(inferno(50)),
                 fontsize=15)
  dev.off()
  cors <- cor(selected[cell_metadata$CellType_reclustering%in%clust, selected_amplicons])
  corrplot::corrplot(cors, diag=FALSE, type='upper')
  
  for_mutual <- selected[cell_metadata$CellType_reclustering%in%clust, selected_amplicons]
  infos <- apply(for_mutual, 2, function(x){
    apply(for_mutual, 2, function(y){
      #sum(x==y)/length(x)
      abs(cor(x, y))
    })
  })
  mean_cors <- apply(infos, 2, mean)
  infos <- infos[order(mean_cors), order(mean_cors)]
  pdf(file.path(plot_path, paste0('correlation_plot_', clust, '.pdf')))
  corrplot::corrplot(infos,
                     tl.col = "black",
                     is.corr = FALSE,
                     type = 'upper',
                     col = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                                                                          "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                                                                          "#4393C3", "#2166AC", "#053061"))(200)))
  dev.off()
  
  pdf(file.path(plot_path, paste0('correlation_heatmap_', clust, '.pdf')), width=3, height=3)
  pheatmap(infos,
           #cutree_rows = 3,
           #cutree_cols = 3,
           clustering_method = 'ward.D2',
           treeheight_row = 0,
           treeheight_col = 0,
           show_rownames=FALSE,
           show_colnames=FALSE,
           border_color=NA,
           legend=FALSE,
           fontsize=8,
           col = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                                  "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                                  "#4393C3", "#2166AC", "#053061"))(200)))
  dev.off()
}
