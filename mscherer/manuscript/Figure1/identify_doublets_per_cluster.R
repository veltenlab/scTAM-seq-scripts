################### identify_doublets_per_cluster.R ################################
#' This file is used to identify further doublets that show substantial deviation from
#' the expected number of features based on a previous clustering of cells.

library(pheatmap)
color_map <- list(CellType=c('naive'='#fcbd7e',
                             'memory1'='#fc6571',
                             'memory2'='#fc3262'),
                  Bulk=c('Naive B-cell high'='#fcbdbd',
                         'Memory B-cell high'='#fc00bd'),
                  CellType_reclustering=c('naive B-cells'='#fcbd7e',
                             'ns-memory B-cells'='#fc6571',
                             'cs-memory B-cells'='#7f47bd'))
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/'

rowinfo <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/cell_metadata.csv',
                    row.names = 1)
selected_amplicons <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/selected_amplicons.csv',
                               row.names = 1)
filtered.counts <- read.table("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/tsv/BCells_Sample7_70_percent_good_performance.barcode.cell.distribution_with_MCL.tsv", row.names = 1, header=T)
amplicon.info <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv", header=T, row.names = 1)

selected_data <- ifelse(filtered.counts[row.names(rowinfo), row.names(selected_amplicons)]>0, 1, 0)

ph <- pheatmap(selected_data, 
               annotation_col = subset(fpr, select = c("Bulk")),
               annotation_row = subset(rowinfo, select = c("CellType", "Nfeatures", "NonCutAmplicons")), 
               clustering_distance_cols = "binary", 
               clustering_distance_rows = "binary", 
               show_rownames = F, 
               show_colnames = F, 
               cutree_rows = 3, 
               clustering_method = "ward.D2",
               color=rev(inferno(50)),
               annotation_colors = color_map,
               fontsize=15)

sel_dat <- subset(rowinfo, subset=CellType=='naive', select='Nfeatures')
hist(unlist(sel_dat))
pot_doublets <- row.names(sel_dat)[sel_dat$Nfeatures>400]
sel_dat <- subset(rowinfo, subset=CellType=='memory1', select='Nfeatures')
hist(unlist(sel_dat))
pot_doublets <- c(row.names(sel_dat)[sel_dat$Nfeatures>380], pot_doublets)
doublets <- rep('Singlet', nrow(rowinfo))
doublets[row.names(rowinfo)%in%pot_doublets] <- 'Doublet'
rowinfo$Doublet <- doublets
png(file.path(plot.path, 'heatmap_complete_doublets_marked.png'),
    width=1000,
    height=1600)
ph <- pheatmap(selected_data, 
               annotation_col = subset(fpr, select = c("Bulk")),
               annotation_row = subset(rowinfo, select = c("CellType", "Nfeatures", "NonCutAmplicons", "Doublet")), 
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
rowinfo[names(re_clustering), 'CellType_reclustering'] <- clust_ass[re_clustering]
png(file.path(plot.path, 'heatmap_complete_reclustering.png'),
    width=1000,
    height=1600)
ph <- pheatmap(selected_data, 
               annotation_col = subset(fpr, select = c("Bulk")),
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
write.csv(rowinfo[rowinfo$Doublet%in%'Singlet', ], file.path(plot.path, 'rowinfo_reclustering_doublet.csv'))
write.csv(pot_doublets, file.path(plot.path, 'reclustering_doublet_names.csv'))