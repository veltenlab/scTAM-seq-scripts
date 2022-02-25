################### F1_D_heatmap.R ################################
#' This file generates the binary DNAm heatmap.

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
uncut <- 'Sample12_70_percent_good_performance'
protein <- 'Sample8'
plot_path <- '~'

selected_amplicons <- read.table(paste0('../../misc/',
                                       uncut, '/tsv/selected_amplicons.tsv'),
                              row.names = 1)
rowinfo <- read.csv(paste0('../../misc/', sample, '/tsv/rowinfo.csv'),
                    row.names = 1)
filtered.counts <- read.table(paste0(sample, "/tsv/", sample, ".barcode.cell.distribution.tsv"),
                                     row.names = 1, header=T)

bulk.methylation <- read.table("../../misc/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt", header=T)
row.names(bulk.methylation) <- bulk.methylation$amplicon
colinfo <- bulk.methylation[row.names(selected_amplicons), c('NBC.mean', 'MBC.mean')]
colinfo$Class <- ifelse(colinfo$NBC.mean>colinfo$MBC.mean, 'Naive B-cell high', 'Memory B-cell high')
selected_data <- ifelse(filtered.counts[row.names(rowinfo), row.names(selected_amplicons)]>0, 1, 0)
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
