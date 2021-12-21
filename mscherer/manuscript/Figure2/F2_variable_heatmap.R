############## F2_variable_heatmap.R ##############
#' With this script, you can plot the heatmap of DNA methylation readouts for the highly variable amplicons.

library(ggplot2)
cut <- 'Sample7_70_percent_good_performance_HP'
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/variable_heatmap/'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=10),
                    axis.text=element_text(color='black',size=8),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')
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
dat_cut <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/", cut, "/tsv/", cut, ".barcode.cell.distribution.tsv"), 
                      sep="\t", 
                      header=T)
cell_metadata <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', 
                            cut, '/tsv/rowinfo.csv'),
                          row.names=1)
amplicon_info <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv",
                          row.names = 1)
amplicon.info <- subset(amplicon.info, subset = Type.of.amplicon=='CpG.B.cell.diff')
dat_cut <- ifelse(dat_cut[row.names(cell_metadata), row.names(amplicon_info)]>0, 1, 0)
for(clust in unique(cell_metadata$CellType_reclustering)){
  cluster_cpgs <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/amplicon_correlation/variable_amplicons_', clust, '.csv'))
  sel_cells <- row.names(cell_metadata)[cell_metadata$CellType_reclustering%in%clust]
  sel_data <- dat_cut[sel_cells, cluster_cpgs$Amplicon]
  png(file.path(plot_path, paste0('heatmap_', clust, '.png')),
      width=1600,
      height=1000)
  ph <- pheatmap(t(sel_data), 
  #               annotation_row = subset(amplicon_info, select='Bulk'),
                 annotation_col = subset(cell_metadata, select = c("Nfeatures")), 
                 clustering_distance_cols = "correlation", 
                 clustering_distance_rows = "correlation", 
                 show_rownames = F, 
                 show_colnames = F, 
                 clustering_method = "ward.D2",
                 color=rev(inferno(50)),
                 annotation_colors = color_map,
                 fontsize=15)
  dev.off()
}