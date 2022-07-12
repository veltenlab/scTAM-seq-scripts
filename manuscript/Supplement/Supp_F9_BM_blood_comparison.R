library(ggplot2)
plot_path <- '~'
blood <- 'GSM5935918_Blood_HhaI'
bm <- 'GSM5935921_BM_HhaI'  
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_line(color='black', size=.25),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')
cluster_blood <- read.csv('../../dropout_modeling', blood, 'corrected_values.csv',
                          row.names=1)
cluster_bm <- read.csv('../../dropout_modeling', bm, 'corrected_values.csv',
                       row.names=1)
joint_names <- intersect(row.names(cluster_blood),
                         row.names(cluster_bm))
cluster_blood <- cluster_blood[joint_names, ]
cluster_bm <- cluster_bm[joint_names, ]
cor_naive <- cor(cluster_blood$Cluster1, cluster_bm$Naive)
cor_cs_memory <- cor(cluster_blood$Cluster2b, cluster_bm$csMBC)
cor_ncs_memory <- cor(cluster_blood$Cluster2a, cluster_bm$ncsMBC)
to_plot <- data.frame(CellType=factor(c(rep('NaiveB', length(joint_names)),
                               rep('CSMemoryB', length(joint_names)),
                               rep('NCSMemoryB', length(joint_names))),levels=c('NaiveB', 'NCSMemoryB', 'CSMemoryB')), 
                      Blood=c(cluster_blood$Cluster1, cluster_blood$Cluster2b, cluster_blood$Cluster2a),
                      BM=c(cluster_bm$Naive, cluster_bm$csMBC, cluster_bm$ncsMBC),
                      Correlation=c(rep(cor_naive, 424),
                                    rep(cor_cs_memory, 424),
                                    rep(cor_ncs_memory, 424)))
color_map <- c('S1cells'='#bdfc7e',
               'S2cells'='#d7ef7e',
               'S3S4cells'='#fcef7e',
               'NaiveB'='#fcbd7e',
               'MemoryB'='#fc6571',
               'NCSMemoryB'='#fc3262',
               'CSMemoryB'='#8e008e')
plot <- ggplot(to_plot, aes(x=Blood, y=BM, color=CellType, group=CellType))+geom_point(size=.1)+geom_smooth(method='lm',se=FALSE, size=.5, color='#808080')+xlim(0, 1)+ylim(0,1)+
  facet_wrap(CellType~.)+
  geom_text(aes(label=paste('r:', format(Correlation, digits=2))), x=0.2, y=0.9, check_overlap=TRUE, color='black', size=1.75, fontface = "bold")+
  plot_theme+
  scale_color_manual(values=color_map)
ggsave(file.path(plot_path, 'memory_naive_B_cells_correlation.pdf'), plot, width=125, height=50, units='mm', dpi=300)
