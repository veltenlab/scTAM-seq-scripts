#### F1_F_G_protein_analysis.R ###################################
#' This file is used to generate plots associated with the protein data of scTAM-seq
#' 

library(ggplot2)
plot_path <- '~'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
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
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_line(color='black', size=.25),
                    strip.background = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')
color_map <- c('naive B-cells'='#fcbd7e',
               'ns-memory B-cells'='#fc3262',
               'cs-memory B-cells'='#8e008e')

dat <- read.csv('../../misc/Sample8_70_percent_good_performance/tsv/rowinfo.csv',
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

thres.cd27 <- 1
thres.cd25 <- 0.5
thres.cd1c <- 1
thres.cd38 <- 1.5
thres.cd71 <- 1.5

is_naive <- sel_dat$CellType%in%'naive B-cells'
p_vals <- apply(sel_dat[, -1], 2, function(x){
  cs <- na.omit(x[is_naive])
  ns <- na.omit(x[!is_naive])
  wilcox.test(cs, ns)$p.value
})
p_vals <- sort(p_vals)

plot <- ggplot(sel_dat, aes(x=CD19, y=CD27, color=CellType))+geom_point(size=.1)+scale_color_manual(values=color_map)+plot_scatter
ggsave(file.path(plot_path, 'antibody_scatterplots_CD19_CD27.pdf'), plot,  width=40, height=40, unit='mm')

q1 <- sel_dat$CD25>thres.cd25&sel_dat$CD27>thres.cd27
q2 <- sel_dat$CD25>thres.cd25&sel_dat$CD27<thres.cd27
q3 <- sel_dat$CD25<thres.cd25&sel_dat$CD27<thres.cd27
q4 <- sel_dat$CD25<thres.cd25&sel_dat$CD27>thres.cd27
frac.q1 <- format(sum(grepl('naive B-cells', sel_dat[q1, 'CellType']))/sum(q1), digit=2)
frac.q2 <- format(sum(grepl('naive B-cells', sel_dat[q2, 'CellType']))/sum(q2), digit=2)
frac.q3 <- format(sum(grepl('naive B-cells', sel_dat[q3, 'CellType']))/sum(q3), digit=2)
frac.q4 <- format(sum(grepl('naive B-cells', sel_dat[q4, 'CellType']))/sum(q4), digit=2)
plot <- ggplot(sel_dat, aes(x=CD27, y=CD25, color=CellType))+geom_point(size=.01)+geom_vline(xintercept = thres.cd27, size=.25)+geom_hline(yintercept = thres.cd25, size=.25)+
  annotate(x=5, y=5, label=frac.q1, geom='text', size=2)+annotate(x=0, y=5, label=frac.q2, geom='text', size=2)+annotate(x=0, y=-1, label=frac.q3, geom='text', size=2)+annotate(x=5, y=-1, label=frac.q4, geom='text', size=2)+
  xlim(-1.5, 6.1)+ylim(-1.5, 6.1)+
    scale_color_manual(values=color_map)+plot_scatter
ggsave(file.path(plot_path, 'antibody_scatterplots_CD27_CD25.pdf'), plot,  width=30, height=30, unit='mm')

to_plot <- reshape2::melt(sel_dat[, !grepl('_raw', colnames(sel_dat))], id='CellType')
colnames(to_plot)[2:3] <- c('Antibody', 'Expression')
plot <- ggplot(to_plot, aes(x=CellType, y=Expression, fill=CellType))+geom_violin()+geom_boxplot(alpha=.25, size=.5, outlier.size = .5)+
  facet_wrap(Antibody~., ncol=8, scale='free_y')+plot_theme+scale_fill_manual(values=color_map)
ggsave(file.path(plot_path, 'antibody_boxplots_clusters_wo_naive_CLR.pdf'), plot)

q1 <- sel_dat$CD38>thres.cd38&sel_dat$CD27>thres.cd27
q2 <- sel_dat$CD38>thres.cd38&sel_dat$CD27<thres.cd27
q3 <- sel_dat$CD38<thres.cd38&sel_dat$CD27<thres.cd27
q4 <- sel_dat$CD38<thres.cd38&sel_dat$CD27>thres.cd27
frac.q1 <- format(sum(grepl('naive B-cells', sel_dat[q1, 'CellType']))/sum(q1), digit=2)
frac.q2 <- format(sum(grepl('naive B-cells', sel_dat[q2, 'CellType']))/sum(q2), digit=2)
frac.q3 <- format(sum(grepl('naive B-cells', sel_dat[q3, 'CellType']))/sum(q3), digit=2)
frac.q4 <- format(sum(grepl('naive B-cells', sel_dat[q4, 'CellType']))/sum(q4), digit=2)
plot <- ggplot(sel_dat, aes(x=CD27, y=CD38, color=CellType))+geom_point(size=.05)+geom_vline(xintercept = thres.cd27)+geom_hline(yintercept = thres.cd38)+
  annotate(x=5, y=5, label=frac.q1, geom='text', size=2)+annotate(x=-1, y=5, label=frac.q2, geom='text', size=2)+annotate(x=-1, y=-1, label=frac.q3, geom='text', size=2)+annotate(x=5, y=-1, label=frac.q4, geom='text', size=2)+
  xlim(-1.5, 6.1)+
  scale_color_manual(values=color_map)+plot_scatter
ggsave(file.path(plot_path, 'antibody_scatterplots_CD27_CD38.pdf'), plot,  width=40, height=40, unit='mm')

q1 <- sel_dat$CD1c>thres.cd1c&sel_dat$CD27>thres.cd27
q2 <- sel_dat$CD1c>thres.cd1c&sel_dat$CD27<thres.cd27
q3 <- sel_dat$CD1c<thres.cd1c&sel_dat$CD27<thres.cd27
q4 <- sel_dat$CD1c<thres.cd1c&sel_dat$CD27>thres.cd27
frac.q1 <- format(sum(grepl('naive B-cells', sel_dat[q1, 'CellType']))/sum(q1), digit=2)
frac.q2 <- format(sum(grepl('naive B-cells', sel_dat[q2, 'CellType']))/sum(q2), digit=2)
frac.q3 <- format(sum(grepl('naive B-cells', sel_dat[q3, 'CellType']))/sum(q3), digit=2)
frac.q4 <- format(sum(grepl('naive B-cells', sel_dat[q4, 'CellType']))/sum(q4), digit=2)
plot <- ggplot(sel_dat, aes(x=CD27, y=CD1c, color=CellType))+geom_point(size=.05)+geom_vline(xintercept = thres.cd27)+geom_hline(yintercept = thres.cd1c)+
  annotate(x=5, y=5, label=frac.q1, geom='text', size=2)+annotate(x=-1, y=5, label=frac.q2, geom='text', size=2)+annotate(x=-1, y=-1, label=frac.q3, geom='text', size=2)+annotate(x=5, y=-1, label=frac.q4, geom='text', size=2)+
  xlim(-1.5, 6.1)+
  scale_color_manual(values=color_map)+plot_scatter
ggsave(file.path(plot_path, 'antibody_scatterplots_CD27_CD1c.pdf'), plot,  width=40, height=40, unit='mm')

q1 <- sel_dat$CD71>thres.cd71&sel_dat$CD27>thres.cd27
q2 <- sel_dat$CD71>thres.cd71&sel_dat$CD27<thres.cd27
q3 <- sel_dat$CD71<thres.cd71&sel_dat$CD27<thres.cd27
q4 <- sel_dat$CD71<thres.cd71&sel_dat$CD27>thres.cd27
frac.q1 <- format(sum(grepl('naive B-cells', sel_dat[q1, 'CellType']))/sum(q1), digit=2)
frac.q2 <- format(sum(grepl('naive B-cells', sel_dat[q2, 'CellType']))/sum(q2), digit=2)
frac.q3 <- format(sum(grepl('naive B-cells', sel_dat[q3, 'CellType']))/sum(q3), digit=2)
frac.q4 <- format(sum(grepl('naive B-cells', sel_dat[q4, 'CellType']))/sum(q4), digit=2)
plot <- ggplot(sel_dat, aes(x=CD27, y=CD71, color=CellType))+geom_point(size=.01)+geom_vline(xintercept = thres.cd27,size=.25)+geom_hline(yintercept = thres.cd71,size=.25)+
  annotate(x=5, y=5, label=frac.q1, geom='text', size=2)+annotate(x=0, y=5, label=frac.q2, geom='text', size=2)+annotate(x=0, y=-1, label=frac.q3, geom='text', size=2)+annotate(x=5, y=-1, label=frac.q4, geom='text', size=2)+
  xlim(-1.5, 6.1)+ylim(-1.5, 6.1)+
  scale_color_manual(values=color_map)+plot_scatter
ggsave(file.path(plot_path, 'antibody_scatterplots_CD27_CD71.pdf'), plot,  width=27, height=27, unit='mm')

is_switched <- sel_dat$CellType%in%'cs-memory B-cells'
sel_dat <- sel_dat[, !grepl('_raw', colnames(sel_dat))]
p_vals <- apply(sel_dat[, -1], 2, function(x){
  cs <- na.omit(x[is_switched])
  ns <- na.omit(x[!is_switched])
  wilcox.test(cs, ns)$p.value
})
p_vals <- sort(p_vals)

sel_dat <- sel_dat[, c('CellType', names(p_vals)[1:6])]
colnames(sel_dat)[-1] <- paste0(colnames(sel_dat)[-1], ', p-value: ', format(p_vals[1:6], digit=3))
to_plot <- reshape2::melt(sel_dat, id='CellType')
colnames(to_plot)[2:3] <- c('Antibody', 'Expression')
plot <- ggplot(to_plot, aes(x=CellType, y=Expression, fill=CellType))+geom_violin()+geom_boxplot(alpha=.25, size=.5, outlier.size = .5)+
  facet_wrap(Antibody~., ncol=3, scale='free_y')+plot_theme+scale_fill_manual(values=color_map)
ggsave(file.path(plot_path, 'antibody_boxplots_clusters_wo_naive_selected_CLR.pdf'), plot,
       width=80, height=40, unit='mm')

sel_dat <- sel_dat[, c('CellType', names(p_vals)[p_vals<0.01])]
colnames(sel_dat)[-1] <- paste0(colnames(sel_dat)[-1], ', p-value: ', format(p_vals[p_vals<0.01], digit=3))
to_plot <- reshape2::melt(sel_dat, id='CellType')
colnames(to_plot)[2:3] <- c('Antibody', 'Expression')
plot <- ggplot(to_plot, aes(x=CellType, y=Expression, fill=CellType))+geom_violin()+geom_boxplot(alpha=.25, size=.5, outlier.size = .5)+
  facet_wrap(Antibody~., ncol=3, scale='free_y')+plot_theme+scale_fill_manual(values=color_map)
ggsave(file.path(plot_path, 'antibody_boxplots_clusters_wo_naive_selected_CLR_significant.pdf'), plot,
       width=80, height=100, unit='mm')

sel_dat <- dat[, !(colnames(dat)%in%c('DoubletDetectionScore',
                                      'DoubletDetectionLabel',
                                      'DoubletClustering',
                                      'Doublet',
                                      'Cluster',
                                      'CellType_broad',
                                      'Nfeatures',
                                      'CD10.1',
                                      'CD117.1'))]
sel_dat <-subset(sel_dat, subset=CellType!='naive B-cells')
sel_dat <-subset(sel_dat, subset=CD25>thres.cd25&CD27>thres.cd27)

q1 <- sel_dat$CD71>thres.cd71&sel_dat$CD1c>thres.cd1c
q2 <- sel_dat$CD71>thres.cd71&sel_dat$CD1c<thres.cd1c
q3 <- sel_dat$CD71<thres.cd71&sel_dat$CD1c<thres.cd1c
q4 <- sel_dat$CD71<thres.cd71&sel_dat$CD1c>thres.cd1c
frac.q1 <- format(sum(grepl('cs-memory B-cells', sel_dat[q1, 'CellType']))/sum(q1), digit=2)
frac.q2 <- format(sum(grepl('cs-memory B-cells', sel_dat[q2, 'CellType']))/sum(q2), digit=2)
frac.q3 <- format(sum(grepl('cs-memory B-cells', sel_dat[q3, 'CellType']))/sum(q3), digit=2)
frac.q4 <- format(sum(grepl('cs-memory B-cells', sel_dat[q4, 'CellType']))/sum(q4), digit=2)
plot <- ggplot(sel_dat, aes(x=CD1c, y=CD71, color=CellType))+geom_point(size=.05)+geom_vline(xintercept = thres.cd1c)+geom_hline(yintercept = thres.cd71)+
  annotate(x=5, y=5, label=frac.q1, geom='text', size=2)+annotate(x=0, y=5, label=frac.q2, geom='text', size=2)+annotate(x=0, y=0, label=frac.q3, geom='text', size=2)+annotate(x=5, y=0, label=frac.q4, geom='text', size=2)+
  xlim(-1.5, 6.1)+ylim(-1.5, 6.1)+
  scale_color_manual(values=color_map)+plot_scatter
ggsave(file.path(plot_path, 'antibody_scatterplots_CD1c_CD71_memory_Bcells.pdf'), plot,  width=27, height=27, unit='mm')

thres.cd27 <- 2
q1 <- sel_dat$CD71>thres.cd71&sel_dat$CD27>thres.cd27
q2 <- sel_dat$CD71>thres.cd71&sel_dat$CD27<thres.cd27
q3 <- sel_dat$CD71<thres.cd71&sel_dat$CD27<thres.cd27
q4 <- sel_dat$CD71<thres.cd71&sel_dat$CD27>thres.cd27
frac.q1 <- format(sum(grepl('cs-memory B-cells', sel_dat[q1, 'CellType']))/sum(q1), digit=2)
frac.q2 <- format(sum(grepl('cs-memory B-cells', sel_dat[q2, 'CellType']))/sum(q2), digit=2)
frac.q3 <- format(sum(grepl('cs-memory B-cells', sel_dat[q3, 'CellType']))/sum(q3), digit=2)
frac.q4 <- format(sum(grepl('cs-memory B-cells', sel_dat[q4, 'CellType']))/sum(q4), digit=2)
plot <- ggplot(sel_dat, aes(x=CD27, y=CD71, color=CellType))+geom_point(size=.05)+geom_vline(xintercept = thres.cd27)+geom_hline(yintercept = thres.cd71)+
  annotate(x=5, y=5, label=frac.q1, geom='text', size=2)+annotate(x=-1, y=5, label=frac.q2, geom='text', size=2)+annotate(x=-1, y=-1, label=frac.q3, geom='text', size=2)+annotate(x=5, y=-1, label=frac.q4, geom='text', size=2)+
  xlim(-1.5, 6.1)+
  scale_color_manual(values=color_map)+plot_scatter
ggsave(file.path(plot_path, 'antibody_scatterplots_CD27_CD71_positive.pdf'), plot,  width=35, height=35, unit='mm')

dat <- read.csv('../../misc/Sample8_70_percent_good_performance/tsv/rowinfo.csv',
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
wilcox.test(sel_dat[sel_dat$CellType%in%'naive B-cells', 'CD27'], sel_dat[!(sel_dat$CellType%in%'naive B-cells'), 'CD27'])
wilcox.test(sel_dat[sel_dat$CellType%in%'naive B-cells', 'CD25'], sel_dat[!(sel_dat$CellType%in%'naive B-cells'), 'CD25'])
wilcox.test(sel_dat[sel_dat$CellType%in%'naive B-cells', 'CD1c'], sel_dat[!(sel_dat$CellType%in%'naive B-cells'), 'CD1c'])
wilcox.test(sel_dat[sel_dat$CellType%in%'naive B-cells', 'CD71'], sel_dat[!(sel_dat$CellType%in%'naive B-cells'), 'CD71'])
wilcox.test(sel_dat[sel_dat$CellType%in%'ns-memory B-cells', 'CD27'], sel_dat[sel_dat$CellType%in%'cs-memory B-cells', 'CD27'])
wilcox.test(sel_dat[sel_dat$CellType%in%'ns-memory B-cells', 'CD25'], sel_dat[sel_dat$CellType%in%'cs-memory B-cells', 'CD25'])
wilcox.test(sel_dat[sel_dat$CellType%in%'ns-memory B-cells', 'CD1c'], sel_dat[sel_dat$CellType%in%'cs-memory B-cells', 'CD1c'])
wilcox.test(sel_dat[sel_dat$CellType%in%'ns-memory B-cells', 'CD71'], sel_dat[sel_dat$CellType%in%'cs-memory B-cells', 'CD71'])
to_plot <- reshape2::melt(sel_dat[, c('CellType', 'CD27', 'CD25', 'CD1c', 'CD71')], id='CellType')
colnames(to_plot)[2:3] <- c('Antibody', 'Expression')
plot <- ggplot(to_plot, aes(x=CellType, y=Expression, fill=CellType))+geom_violin(size=.25)+geom_boxplot(alpha=.25, size=.25, outlier.size = .25)+
  facet_wrap(Antibody~., ncol=4, scale='free_y')+plot_theme+scale_fill_manual(values=color_map)+ylab('Surface Expression')
ggsave(file.path(plot_path, 'antibody_boxplots_clusters_selected.pdf'), plot, width=120, height=40, unit='mm')
