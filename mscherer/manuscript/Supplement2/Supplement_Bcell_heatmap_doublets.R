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
sample <- 'Sample8_70_percent_good_performance'
filtered.counts <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/", sample, "/tsv/", sample, ".barcode.cell.distribution.tsv"), 
                              sep="\t", 
                              header=T)
amplicon.info <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv", header=T, row.names = 1)

bulk.methylation <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt", header=T)
rowinfo <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/all_doublets.tsv'),
                    row.names = 1)
colinfo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/colinfo.csv'
                           ), row.names=1)
png(file.path(plot_path, 'heatmap_doublet.png'),
    width=1000,
    height=1600)
ph <- pheatmap(ifelse(filtered.counts[row.names(rowinfo), row.names(colinfo)]>0, 1, 0), 
               #annotation_col = colinfo,
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
