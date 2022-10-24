#### F2_D_protein_analysis.R ###################################
#' This file is used to generate plots associated with the protein data of scTAM-seq
#' 

library(ggplot2)
plot_path <- '~'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=7),
                    axis.text=element_text(color='black',size=6),
                    axis.ticks=element_line(color='black', size=.25),
                    strip.background = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    #axis.text.x=element_text(angle=45, hjust=1, vjust = 1),
                    axis.text.x=element_blank(),
                    axis.title.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    legend.position='none',
                    panel.spacing = unit(.1, "lines"))
plot_scatter <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=7),
                    axis.text=element_text(color='black',size=6),
                    axis.ticks=element_line(color='black', size=.25),
                    strip.background = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')
color_map <- c('naive B-cells'='#fcbd7e',
               'ns-memory B-cells'='#fc3262',
               'cs-memory B-cells'='#8e008e')

dat <- read.csv('../../GSM5935918_Blood_HhaI/tsv/rowinfo_NCS.csv',
                    row.names = 1)
dat$CellType <- factor(as.character(dat$CellType), levels=c('naive B-cells',
                                                            'ns-memory B-cells',
                                                            'cs-memory B-cells'))
sel_dat <- dat[, !(colnames(dat)%in%c('DoubletDetectionScore',
                                              'DoubletDetectionLabel',
                                              'DoubletClustering',
                                              'Doublet',
                                              'Cluster',
                                              'CellType_broad',
                                              'Nfeatures',
                                              'CD10.1',
                                              'CD117.1'))]
to_plot <- reshape2::melt(sel_dat[, !grepl('_raw', colnames(sel_dat))], id='CellType')
colnames(to_plot)[2:3] <- c('Antibody', 'Expression')
plot <- ggplot(to_plot, aes(x=CellType, y=Expression, fill=CellType))+geom_violin(size=.25)+geom_boxplot(alpha=.25, size=.25, outlier.size = .25)+
  facet_wrap(Antibody~., ncol=12, scale='free_y')+plot_theme+scale_fill_manual(values=color_map)
ggsave(file.path(plot_path, 'antibody_boxplots_clusters_CLR.pdf'), plot, width=175, height=65, unit='mm')

dat$ncBCellClustering_detailed <- dat$ncBCellClustering
dat[dat$ncBCellClustering_detailed%in%'Other', 'ncBCellClustering_detailed'] <- as.character(dat$CellType[dat$ncBCellClustering_detailed%in%'Other'])
dat$ncBCellClustering_detailed <- factor(c('1'='Cluster 1',
                                    '2_1'='Cluster 2_1',
                                    '2_2'='Cluster 2_2',
                                    '4'='Cluster 4',
                                    '5'='Cluster 5',
                                    'naive B-cells'='naive B-cells',
                                    'cs-memory B-cells'='cs-memory B-cells')[as.character(dat$ncBCellClustering_detailed)],
levels=c('naive B-cells', 'Cluster 1', 'Cluster 2_2', 'Cluster 2_1', 'Cluster 4', 'Cluster 5', 'cs-memory B-cells'))
to_plot <- dat[, c('ncBCellClustering_detailed',
                   'CD27',
                   'CD11c')]
to_plot <- reshape2::melt(to_plot, id='ncBCellClustering_detailed')
colnames(to_plot) <- c('Cluster', 'Antibody', 'Expression')
color_map2 <- c('naive B-cells'='#fcbd7e',
                'Cluster 1'='#4a0f1d',
                'Cluster 2_2'='#843c15',
                'Cluster 2_1'='#c45a22',
                'Cluster 4'='#f2615f',
                'Cluster 5'='#fbe5d9',
               'cs-memory B-cells'='#8e008e')
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
plot <- ggplot(to_plot, aes(x=Cluster, y=Expression, fill=Cluster))+geom_violin(size=.25)+geom_boxplot(alpha=.25, size=.25, outlier.size = .25)+
  facet_wrap(Antibody~., ncol=2, scale='free_y')+plot_theme2+scale_fill_manual(values=color_map2)
ggsave(file.path(plot_path, 'antibody_boxplots_clusters_memory.pdf'), plot, width=45, height=25, unit='mm')

to_plot <- dat[, c('CD27',
                   'CD11c')]
mod <- data.frame(to_plot, Cluster=factor(as.character(dat$ncBCellClustering_detailed),
                                  levels=c('naive B-cells', 'Cluster 1', 'Cluster 2_2', 'Cluster 2_1', 'Cluster 4', 'Cluster 5', 'cs-memory B-cells')))
lm.mod <- lm(CD27~Cluster, data=mod)
kruskal.test(mod$CD27, g=mod$Cluster)
lm.mod <- lm(CD11c~Cluster, data=mod)
kruskal.test(mod$CD11c, g=mod$Cluster)
ncBCellClustering_detailed <- as.character(dat$ncBCellClustering_detailed)
to_plot <- aggregate(to_plot, by=list(ncBCellClustering_detailed), function(x){c(mean(x), sd(x)/sqrt(length(x)))})
to_plot <- data.frame(Cluster=rep(to_plot$Group.1, 2),
                      Antibody=c(rep('CD27', 7), rep('CD11', 7)),
                      Mean=c(to_plot$CD27[, 1], to_plot$CD11[, 1]),
                      SD=c(to_plot$CD27[, 2], to_plot$CD11[, 2]))
to_plot$Cluster <- factor(to_plot$Cluster,
                          levels=c('naive B-cells', 'Cluster 1', 'Cluster 2_2', 'Cluster 2_1', 'Cluster 4', 'Cluster 5', 'cs-memory B-cells'))

plot <- ggplot(to_plot, aes(x=Cluster, y=Mean, fill=Cluster))+geom_bar(stat='identity')+geom_errorbar(aes(ymax=Mean+2*SD, ymin=Mean-2*SD), size=.25, width=.5)+
  facet_wrap(Antibody~., ncol=2, scale='free_y')+plot_theme2+scale_fill_manual(values=color_map2)+scale_color_manual(values=color_map2)
ggsave(file.path(plot_path, 'antibody_barlots_clusters_memory.pdf'), plot, width=80, height=25, unit='mm')
