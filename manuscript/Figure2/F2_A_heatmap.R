library(ComplexHeatmap)
library(viridis)
library(ggplot2)
color_map <- list(cluster=c('1'='gray10',
                            '2'='gray23',
                            '3'='gray41',
                            '4'='gray59',
                            '5'='gray77',
                            '6'='gray93'),
                  ncBCellClustering=c('naive B-cells'='#fcbd7e',
                                      '1'='#4a0f1d',
                                      '2_2'='#843c15',
                                      '2_1'='#c45a22',
                                      '4'='#f2615f',
                                      '5'='#fbe5d9',
                                      'cs-memory B-cells'='#8e008e'),
                  Class=c('Naive B-cell high'='#fcbdbd',
                          'Memory B-cell high'='#fc00bd'),
                  ncBCellClustering=c('ns-memory B-cells'='#fc3262',
                             'cs-memory B-cells'='#8e008e'),
                  MBC.mean=rev(c('#43a2ca',
                                 '#7bccc4',
                                 '#bae4bc',
                                 '#f0f9e8')),
                  NBC.mean=rev(c('#8856a7',
                                 '#8c96c6',
                                 '#b3cde3',
                                 '#edf8fb')))
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(angle=90, hjust=1, color='black',size=5),
                    axis.ticks=element_line(color='black', size=.25),
                    strip.background = element_blank(),
                    #strip.text.x = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')
sample <- 'GSM5935918_Blood_HhaI'
uncut <- 'GSM5935923_BM_undigested'
plot_path <- '~'

selected_amplicons <- read.table(paste0('../../misc/',
                                        uncut, '/tsv/selected_amplicons.tsv'),
                                 row.names = 1)
rowinfo <- read.csv(paste0('../../misc/', sample, '/tsv/rowinfo.csv'),
                    row.names = 1)
filtered.counts <- read.table(paste0("../../data/", sample, ".tsv.gz"),
                              row.names = 1, header=T)

bulk.methylation <- read.table("../../misc/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt",
                               header=T)
row.names(bulk.methylation) <- bulk.methylation$amplicon
colinfo <- bulk.methylation[row.names(selected_amplicons), c('NBC.mean', 'MBC.mean')]
colinfo$Class <- ifelse(colinfo$NBC.mean>colinfo$MBC.mean, 'Naive B-cell high', 'Memory B-cell high')
selected_data <- ifelse(filtered.counts[row.names(rowinfo), row.names(selected_amplicons)]>0, 1, 0)

selected_data <- selected_data[row.names(subset(rowinfo, subset = CellType=='ns-memory B-cells')), ]

clustering <- cutree(hclust(dist(selected_data, 'binary'), method='ward.D2'), 5)
clustering_names <- c('1'='5',
                      '2'='1',
                      '3'='4',
                      '4'='2_1',
                      '5'='2_2')[clustering]
rowinfo$ncBCellClustering <- 'Other'
rowinfo[names(clustering), 'ncBCellClustering'] <- clustering_names
sel_cells <- names(clustering[clustering==2])
counts <- filtered.counts[row.names(rowinfo), row.names(selected_amplicons)]
counts_c1 <- counts[sel_cells, ]
counts_c2 <- counts[!(row.names(counts)%in%sel_cells), ]
p_vals <- c()
for(i in 1:ncol(counts_c1)){
  log2_fc <- log2(mean(counts_c1[, i])/mean(counts_c2[, i]))
  p_val <- wilcox.test(counts_c1[, i], counts_c2[, i])$p.value
  p_vals <- rbind(p_vals, data.frame(p_val, log2_fc))
}
p_vals$p.adj <- p.adjust(p_vals$p_val)
row.names(p_vals) <- colnames(counts_c1)
colinfo$logP <- -log10(p_vals$p_val)
colinfo$selected <- ifelse(p_vals$p.adj<0.01, 'Yes', 'No')
rowinfo <- subset(rowinfo, subset = CellType=='ns-memory B-cells')
clust_ampls <- as.factor(cutree(hclust(dist(t(selected_data), 'binary'), method='ward.D2'), 6))
colinfo$cluster <- clust_ampls
color_map <- list(c)
png(file.path(plot_path, 'heatmap_ns_B_cells.png'),
    width=3000,
    height=4800,
    res=300)
pheatmap(selected_data, 
         annotation_col = subset(colinfo, select = c("NBC.mean",
                                                     "MBC.mean",
                                                     'cluster')),
         annotation_row = subset(rowinfo, select = c("ncBCellClustering")),
         clustering_distance_cols = "binary", 
         clustering_distance_rows = "binary", 
         show_rownames = F, 
         show_colnames = F, 
         cutree_rows = 5, 
         cutree_cols = 6, 
         clustering_method = "ward.D2",
         color=rev(inferno(50)),
         annotation_colors = color_map,
         fontsize=15,
         labels_row = FALSE,
         labels_col = FALSE,
         raster_by_magick = TRUE,
         raster_quality=10)
dev.off()
write.csv(colinfo, '../../misc/colinfo_diff_cells.csv')

clust.means <- aggregate(selected_data, by=list(rowinfo$ncBCellClustering), mean)
row.names(clust.means) <- clust.means$Group.1
clust.means <- clust.means[, -1]
clust.sds <- aggregate(t(clust.means), by=list(colinfo$cluster), function(x){
  sd(x)/sqrt(length(x))
})
clust.means <- aggregate(t(clust.means), by=list(colinfo$cluster), mean)
to_plot <- data.frame(reshape2::melt(clust.means, id='Group.1'),
                      SE=reshape2::melt(clust.sds, id='Group.1')$value)
colnames(to_plot) <- c('CpGCluster', 'CellCluster', 'Methylation', 'SE')
plot_theme2 <- theme(panel.background = element_rect(color='black',fill='white'),
                     panel.grid=element_blank(),
                     text=element_text(color='black',size=6),
                     axis.text=element_text(color='black',size=5),
                     axis.ticks=element_line(color='black', size=.25),
                     strip.background = element_blank(),
                     legend.key=element_rect(color=NA, fill=NA),
                     axis.text.x=element_blank(),
                     axis.title.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     panel.spacing = unit(.1, "lines"),
                     legend.position = 'none')
to_plot$'CellCluster' <- factor(to_plot$'CellCluster',
                                 levels = c('1', '2_2', '2_1', '4', '5'))
plot <- ggplot(to_plot, aes_string(x='CellCluster', y='Methylation', fill='CellCluster'))+geom_bar(stat='identity')+geom_errorbar(aes(ymax=Methylation+2*SE, ymin=Methylation-2*SE), size=.25, width=.5)+
  facet_wrap(Antibody~., ncol=2, scale='free_y')+plot_theme2+scale_fill_manual(values=color_map$ncBCellClustering)+scale_color_manual(values=color_map$ncBCellClustering)+facet_wrap(CpGCluster~., ncol=3)
ggsave(file.path(plot_path, 'cluster_methylation.pdf'), plot, width=75, height=50, unit='mm')
