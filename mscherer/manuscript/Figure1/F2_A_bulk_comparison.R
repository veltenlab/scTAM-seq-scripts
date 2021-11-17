################### F1_E_bulk_comparison.R ################### 
#'This file compares the pseudo-bulk obtained from averaging the methylation calls
#'in the clusters defined in the heatmap with the bulk values measured in
#'doi.org/10.1038/ng.3291


library(ggplot2)
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=10),
                    axis.text=element_text(color='black',size=8),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    #strip.text.x = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5),
               legend.position='none')
color_map <- c('Cluster1'='#fcbd7e',
               'Cluster2a'='#fc3262',
               'Cluster2b'='#bd0071',
               'Cluster2c'='#8e008e')

means_clusters <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/dropou_modeling/all_corrected_stan_clusters.csv', row.names=1)
more_info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
row.names(more_info) <- more_info$amplicon
joint_names <- intersect(row.names(more_info), row.names(means_clusters))
means_clusters <- means_clusters[joint_names, ]
row.names(means_clusters) <- more_info[row.names(means_clusters), 'background.cpgs']

load('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/meth.data.numeric.Rdata')
cell_assignment <- read.table('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/MBC_assignment.txt')
meth.data.numeric <- meth.data.numeric[,as.character(cell_assignment$V5)]
means_csMBC <- apply(meth.data.numeric[row.names(means_clusters), cell_assignment$V2=='csMBC'], 1, mean)
means_ncsMBC <- apply(meth.data.numeric[row.names(means_clusters), cell_assignment$V2=='ncsMBC'], 1, mean)
to_plot <- data.frame(means_clusters, csMBC=means_csMBC, ncsMBC=means_ncsMBC)
to_plot$Amplicon <- row.names(to_plot)
to_plot <- reshape2::melt(to_plot, id=c('Amplicon', 'csMBC', 'ncsMBC'))
colnames(to_plot)[4:5] <- c('Cluster', 'Methylation')
to_plot <- reshape2::melt(to_plot, id=c('Amplicon', 'Cluster', 'Methylation'))
colnames(to_plot)[4:5] <- c('Bulk', 'BulkMethylation')
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
to_plot$BulkMethylation <- as.numeric(to_plot$BulkMethylation)
cors <- t(as.data.frame(apply(to_plot, 1, function(x){
  ct_sc <- x['Cluster']
  ct_bulk <- x['Bulk']
  meth <- subset(to_plot, Cluster==ct_sc & Bulk==ct_bulk, select=c('Methylation', 'BulkMethylation'))
  c(cor(meth$Methylation, meth$BulkMethylation),
    cor.test(meth$Methylation, meth$BulkMethylation)$p.value)
})))
to_plot$Correlation <- cors[,1]
to_plot$CorrelationPVal <- cors[,2]
plot <- ggplot(to_plot, aes_string(x='BulkMethylation', y='Methylation', color='Cluster'))+geom_point()+geom_smooth(method='lm',se=FALSE)+xlim(0, 1)+ylim(0,1)+
  facet_wrap(Bulk~Cluster, nrow=4)+
  geom_text(aes(label=paste('Pearsons R²:', format(Correlation, digits=2))), x=0.35, y=0.9, check_overlap=TRUE, color='black')+
  plot_theme+
  scale_color_manual(values=color_map)+ylab('Pseudo-bulk methylation')+xlab('Bulk methylation')
ggsave(file.path(plot_path, 'F2_A_bulk_vs_singlecell.pdf'), plot, width=250, height=120, units='mm')
