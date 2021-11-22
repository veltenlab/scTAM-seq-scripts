################### F2_explore_imprinted_behavior.R ################### 
#' This files explores the data for indications of allele-specific methylation and for cellular subpopulations.

library(ggplot2)
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=10),
                    axis.text=element_text(color='black',size=8),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5),
                    legend.position='none')

means_clusters <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/dropou_modeling/all_corrected_stan_clusters.csv', row.names=1)
dat_cut <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/", cut, "/tsv/", cut, ".barcode.cell.distribution_with_MCL.tsv"), 
                      sep="\t", 
                      header=T)
doublet <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/reclustering_doublet_names.csv')
double2 <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/tsv/doublet_scores_DoubletDetection.csv')
doublets <- c(doublet$x, double2$Barcode[double2$DoubletDetectionLabel==1], row.names(means_clusters))
dat_cut <- dat_cut[!(row.names(dat_cut)%in%doublets), ]
rowinfo <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/rowinfo_reclustering_doublet.csv',
                    row.names = 1)
dat_cut <- dat_cut[row.names(rowinfo), row.names(means_clusters)]
meth_clusters <- aggregate(dat_cut, by=list(rowinfo$CellType_reclustering), function(x){
  sum(x>0)/length(x)
})
row.names(meth_clusters) <- meth_clusters$Group.1
meth_clusters <- as.data.frame(t(meth_clusters[, -1]))

dat_uncut <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/", uncut, "/tsv/", uncut, ".barcode.cell.distribution.tsv"), 
                        sep="\t", 
                        header=T)
doublet <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample6_70_percent_good_performance/tsv/doublet_scores_DoubletDetection.csv')
doublets <- c(doublet$Barcode[doublet$DoubletDetectionLabel==1])
dat_uncut <- dat_uncut[!(row.names(dat_uncut)%in%doublets), row.names(means_clusters)]
meth_rates <- apply(ifelse(dat_uncut>0, 1, 0), 2, function(x){
  sum(x>0)/length(x)
})
means_clusters$Uncut <- meth_rates
meth_clusters$Uncut <- meth_rates
to_plot <- reshape2::melt(meth_clusters, id='Uncut')
colnames(to_plot)[2:3] <- c('Cluster', 'Methylation')
distance.from.square.fit <- apply(to_plot, 1, function(x){
  if(x[3]<0.25|x[3]>0.75) return(1)
  abs(sqrt(as.numeric(x[3]))-as.numeric(x[1]))
})
to_plot$Distance <- distance.from.square.fit
color_map <- c('Cluster1'='#fcbd7e',
               'Cluster2a'='#fc3262',
               'Cluster2b'='#bd0071',
               'Cluster2c'='#8e008e')
square <- data.frame(Uncut=sqrt(seq(0,1,by=0.001)), Cut=seq(0,1,by=0.001))
plot <- ggplot()+geom_point(data=square, aes(x=Uncut, y=Cut))+geom_point(data=to_plot, aes(x=Uncut, y=Methylation, color=Cluster, alpha=1-Distance))+geom_abline(slope=1, intercept=0)+
 facet_wrap(Cluster~.)+scale_color_manual(values=color_map)+plot_theme+
  labs(x='Fraction of cells with reads in undigested sample', y='Fraction of cells with reads in digested sample')
ggsave('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/F2_E_imprinted_in_amplicons.pdf', plot,
       width=140,
       height=90,
       units='mm')

mean_distances <- apply(means_clusters, 1, function(x){
  if(any(x[5]<0.25|x[-5]>0.75)) return(1)
  mean(abs(sqrt(as.numeric(x[5]))-as.numeric(x[-5])))
})
