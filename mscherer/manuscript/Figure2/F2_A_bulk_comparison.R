################### F1_E_bulk_comparison.R ################### 
#'This file compares the pseudo-bulk obtained from averaging the methylation calls
#'in the clusters defined in the heatmap with the bulk values measured in
#'doi.org/10.1038/ng.3291


library(ggplot2)
plot_path <- '/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure2/Sample8/'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_line(color='black', size=.25),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5),
               legend.position='none')

means_clusters <- read.csv('/users/lvelten/project/Methylome/analysis/dropou_modeling/re_clustering/Sample8/clusters/all_amplicons_clusters_batch_corrected.csv', row.names=1)
more_info <- read.table('/users/lvelten/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
row.names(more_info) <- more_info$amplicon
joint_names <- intersect(row.names(more_info), row.names(means_clusters))
means_clusters <- means_clusters[joint_names, ]
row.names(means_clusters) <- more_info[row.names(means_clusters), 'background.cpgs']

load('/users/lvelten/project/Methylome/data/external/BLUEPRINT/Renee/meth.data.numeric.Rdata')
cell_assignment <- read.table('/users/lvelten/project/Methylome/data/external/BLUEPRINT/Renee/MBC_assignment.txt')
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
to_plot <- na.omit(to_plot)
cors <- t(as.data.frame(apply(to_plot, 1, function(x){
  ct_sc <- x['Cluster']
  ct_bulk <- x['Bulk']
  meth <- subset(to_plot, Cluster==ct_sc & Bulk==ct_bulk, select=c('Methylation', 'BulkMethylation'))
  c(cor(meth$Methylation, meth$BulkMethylation),
    cor.test(meth$Methylation, meth$BulkMethylation)$p.value)
})))
to_plot$Correlation <- cors[,1]
to_plot$CorrelationPVal <- cors[,2]
color_map <- c('Cluster1'='#fcbd7e',
               'Cluster2a'='#fc3262',
#               'Cluster2b'='#bd0071',
               'Cluster2b'='#8e008e')
plot <- ggplot(to_plot, aes_string(x='BulkMethylation', y='Methylation', color='Cluster'))+geom_point(size=.5)+geom_smooth(method='lm',se=FALSE)+xlim(0, 1)+ylim(0,1)+
  facet_wrap(Cluster~Bulk, nrow=4)+
  geom_text(aes(label=paste('Pearsons R²:', format(Correlation, digits=2))), x=0.3, y=0.9, check_overlap=TRUE, color='black', size=2.5)+
  plot_theme+
  scale_color_manual(values=color_map)+ylab('Pseudo-bulk methylation')+xlab('Bulk methylation')
ggsave(file.path(plot_path, 'F2_A_bulk_vs_singlecell.pdf'), plot, width=90, height=63, units='mm')

color_map <- c('1'='#fcbd7e',
               '2'='#fc3262',
               '3'='#bd0071',
               '4'='#8e008e',
               '5'='#e78ac3',
               '6'=NA)
to_plot <- data.frame(means_clusters, csMBC=means_csMBC, ncsMBC=means_ncsMBC)
to_plot$Amplicon <- row.names(to_plot)
to_plot <- reshape2::melt(to_plot, id=c('Amplicon', 'csMBC', 'ncsMBC'))
colnames(to_plot)[4:5] <- c('Cluster', 'Methylation')
to_plot <- reshape2::melt(to_plot, id=c('Amplicon', 'Cluster', 'Methylation'))
colnames(to_plot)[4:5] <- c('Bulk', 'BulkMethylation')
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
to_plot$BulkMethylation <- as.numeric(to_plot$BulkMethylation)
to_plot$Type <- 'BCellDifferentiation'
imprinted_bulk <- read.csv('/users/lvelten/project/Methylome/infos/BCells/imprinted_amplicons_bulk.csv', row.names=1)
imprinted_clusters <- read.csv('/users/lvelten/project/Methylome/analysis/dropou_modeling/corrected_stan_clusters_imprinted.csv', row.names=1)
joint_names <- intersect(row.names(imprinted_bulk), colnames(imprinted_clusters))
imprinted_bulk <- imprinted_bulk[joint_names, ]
imprinted_clusters <- imprinted_clusters[, joint_names]
to_add <- data.frame(csMBC=imprinted_bulk$ImprintedBulk, ncsMBC=imprinted_bulk$ImprintedBulk, t(imprinted_clusters), Amplicon=joint_names)
to_add <- reshape2::melt(to_add, id=c('Amplicon', 'csMBC', 'ncsMBC'))
colnames(to_add)[4:5] <- c('Cluster', 'Methylation')
to_add <- reshape2::melt(to_add, id=c('Amplicon', 'Cluster', 'Methylation'))
colnames(to_add)[4:5] <- c('Bulk', 'BulkMethylation')
to_add$Correlation <- NA
to_add$CorrelationPVal <- NA
to_add$Type <- 'Imprinted'
to_plot <- as.data.frame(rbind(to_plot, to_add))
to_plot <- to_plot[(to_plot$Methylation>0.25)&(to_plot$Methylation<0.75), ]
plot <- ggplot(to_plot, aes_string(x='BulkMethylation', y='Methylation'))+geom_point(aes(color=ifelse(Type=='Imprinted', '5', Cluster)))+geom_smooth(aes(color=ifelse(Type=='Imprinted', '6', Cluster)), method='lm',se=FALSE)+xlim(0, 1)+ylim(0,1)+
  facet_wrap(Bulk~Cluster, nrow=4)+
  geom_text(aes(label=paste('Pearsons R²:', format(Correlation, digits=2))), x=0.25, y=0.9, check_overlap=TRUE, color='black', size=.5)+
  plot_theme+
  scale_color_manual(values=color_map)+ylab('Pseudo-bulk methylation')+xlab('Bulk methylation')
ggsave(file.path(plot_path, 'F2_A_bulk_vs_singlecell_imprinted.pdf'), plot, width=90, height=63, units='mm')
