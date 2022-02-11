################### Figure_2_heatmap_ncs_memory_B_cells.R ################################
#' With this script, we compare the two subcluster that we identify in ncs-memory B-cells

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
                             'cs-memory B-cells'='#fc3262'))
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_line(color='black', size=.25),
                    strip.background = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    #axis.text.x=element_text(angle=45, hjust=1, vjust = 1),
                    #axis.text.x=element_blank(),
                    #axis.title.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    legend.position='none',
                    panel.spacing = unit(.1, "lines"))
sample <- 'Sample8_70_percent_good_performance'
uncut <- 'Sample6_70_percent_good_performance'
protein <- 'Sample8'
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/Sample8/'

selected_amplicons <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv',
                                 row.names = 1)
selected_amplicons <- subset(selected_amplicons, subset = Type.of.amplicon=="CpG.B.cell.diff")
rowinfo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/rowinfo.csv'),
                    row.names = 1)
filtered.counts <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/", sample, "/tsv/", sample, ".barcode.cell.distribution.tsv"),
                              row.names = 1, header=T)

bulk.methylation <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt", header=T)
row.names(bulk.methylation) <- bulk.methylation$amplicon
colinfo <- bulk.methylation[row.names(selected_amplicons), c('NBC.mean', 'MBC.mean')]
colinfo$Class <- ifelse(colinfo$NBC.mean>colinfo$MBC.mean, 'Naive B-cell high', 'Memory B-cell high')
selected_data <- ifelse(filtered.counts[row.names(rowinfo), row.names(selected_amplicons)]>0, 1, 0)
selected_data <- selected_data[which(rowinfo$Doublet==0), ]
rowinfo <- rowinfo[which(rowinfo$Doublet==0), ]

selected_data <- selected_data[which(rowinfo$Cluster=='Cluster2a'), ]
rowinfo <- rowinfo[which(rowinfo$Cluster=='Cluster2a'), ]

load('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/meth.data.numeric.Rdata')
cell_assignment <- read.table('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/MBC_assignment.txt')
meth.data.numeric <- meth.data.numeric[,as.character(cell_assignment$V5)]
more_info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
row.names(more_info) <- more_info$amplicon
cps <- as.character(more_info[row.names(selected_amplicons), 'background.cpgs'])
mean_classes <- aggregate(t(meth.data.numeric[cps, ]), by=list(cell_assignment$V2), mean)
row.names(mean_classes) <- mean_classes$Group.1
mean_classes <- mean_classes[ ,-1]
colnames(mean_classes) <- row.names(selected_amplicons)
colinfo <- data.frame(colinfo, t(mean_classes))

png(file.path(plot_path, 'heatmap_complete_ncs_memory_Bcells.png'),
    width=1000,
    height=1600)
ph <- pheatmap(selected_data, 
               annotation_col = subset(colinfo, select = c("csMBC",
                                                           "ncsMBC")),
               annotation_row = subset(rowinfo, select = c("CellType",
                                                           "Nfeatures",
                                                           "CD27")),
               clustering_distance_cols = "binary", 
               clustering_distance_rows = "binary", 
               show_rownames = F, 
               show_colnames = F, 
               cutree_rows = 2, 
               clustering_method = "ward.D2",
               color=rev(inferno(50)),
               annotation_colors = color_map,
               fontsize=15)
dev.off()

clust <- hclust(dist(selected_data, 'binary'), 'ward.D2')
clust <- cutree(clust, 2)
rowinfo$ClustNCS <- as.factor(clust)

plot <- ggplot(rowinfo, aes(x=ClustNCS, y=CD27))+geom_boxplot()+plot_theme
ggsave(file.path(plot_path, 'CD27_boxplot_clusters.pdf'), plot, width=100, height=100, unit='mm')

to_plot <- rowinfo[, c(9:56, 103)]
to_plot <- reshape2::melt(to_plot, id='ClustNCS')
colnames(to_plot)[2:3] <- c('Antibody', 'Expression')
plot <- ggplot(to_plot, aes(x=ClustNCS, y=Expression))+geom_violin()+geom_boxplot(alpha=.25, size=.5, outlier.size = .5)+
  facet_wrap(Antibody~., ncol=12, scale='free_y')+plot_theme
ggsave(file.path(plot_path, 'antibody_boxplot_clusters.pdf'), plot)
