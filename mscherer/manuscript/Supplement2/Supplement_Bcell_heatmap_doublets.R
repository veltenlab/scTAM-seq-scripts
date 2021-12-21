############################### Supplement_Bcell_heatmap_doublets.R ######################
#' This file generates the heatmap after filtering for amplicons that work well in the 
#' undigested Sample. 

.libPaths(c(.libPaths(), '/users/mscherer/conda/envs/rnbeads/lib/R/library/'))
library(Seurat)
library(ggplot2)
library(plyr)
library(reshape2)
library(pheatmap)
library(viridis)
library(grid)
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Supplement/'
color_map <- list(CellType=c('naive'='#fcbd7e',
                 'memory1'='#fc6571',
                 'memory2'='#fc3262'),
               Bulk=c('Naive B-cell high'='#fcbd7e',
               'Memory B-cell high'='#fc4762'))
sample <- 'Sample7_70_percent_good_performance'
filtered.counts <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/", sample, "/tsv/", cut, ".barcode.cell.distribution.tsv"), 
                              sep="\t", 
                              header=T)
amplicon.info <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv", header=T, row.names = 1)

bulk.methylation <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt", header=T)
rowinfo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/doublet_scores_DoubletDetection.csv'),
                    row.names = 1)
colinfo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/colinfo.csv'))
png(file.path(plot_path, 'heatmap_doublet.png'),
    width=1000,
    height=1600)
ph <- pheatmap(selected_data, 
               annotation_col = colinfo,
               annotation_row = subset(rowinfo, select = c("DoubletDetectionLabel", "DoubletDetectionScore", "Doublet")), 
               clustering_distance_cols = "binary", 
               clustering_distance_rows = "binary", 
               show_rownames = F, 
               show_colnames = F, 
               clustering_method = "ward.D2",
               color=rev(inferno(50)),
               annotation_colors = color_map,
               fontsize=15)
dev.off()
