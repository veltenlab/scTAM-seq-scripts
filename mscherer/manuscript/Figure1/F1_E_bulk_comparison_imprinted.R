################### F1_E_bulk_comparison.R ################### 
#'This file compares the pseudo-bulk obtained from averaging the methylation calls
#'in the clusters defined in the heatmap with the bulk values measured in
#'doi.org/10.1038/ng.3291


library(ggplot2)
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1//'
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
color_map <- c('Naive'='#fcbd7e',
               'Memory'='#fc3262')

means_clusters <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/dropou_modeling/all_corrected_stan_memory_naive.csv', row.names=1)
bulk_table <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
means_bulk <- bulk_table[, c('NBC.mean', 'MBC.mean')]
colnames(means_bulk) <- c('naiveB', 'memoryB')
row.names(means_bulk) <- bulk_table$amplicon
means_bulk <- as.data.frame(t(means_bulk))
means_bulk$'Group.1' <- row.names(means_bulk)
to_plot <- as.data.frame(cbind(means_clusters, t(means_bulk)[row.names(means_clusters)[row.names(means_clusters)%in%colnames(means_bulk)], ]))
to_plot$Amplicon <- row.names(to_plot)
to_plot <- reshape2::melt(to_plot, id=c('Amplicon', 'naiveB', 'memoryB'))
colnames(to_plot)[4:5] <- c('CellType', 'Methylation')
to_plot <- reshape2::melt(to_plot, id=c('Amplicon', 'CellType', 'Methylation'))
colnames(to_plot)[4:5] <- c('Bulk', 'BulkMethylation')
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
to_plot$BulkMethylation <- as.numeric(to_plot$BulkMethylation)
cors <- t(as.data.frame(apply(to_plot, 1, function(x){
  ct_sc <- x['CellType']
  ct_bulk <- x['Bulk']
  meth <- subset(to_plot, CellType==ct_sc & Bulk==ct_bulk, select=c('Methylation', 'BulkMethylation'))
  c(cor(meth$Methylation, meth$BulkMethylation),
    cor.test(meth$Methylation, meth$BulkMethylation)$p.value)
})))
to_plot$Correlation <- cors[,1]
to_plot$CorrelationPVal <- cors[,2]
to_plot$Type <- 'BCellDifferentiation'
imprinted_bulk <- read.csv('/users/mscherer/cluster/project/Methylome/infos/BCells/imprinted_amplicons_bulk.csv', row.names=1)
imprinted_clusters <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/dropou_modeling/corrected_stan_clusters_imprinted.csv', row.names=1)
imprinted_clusters <- as.data.frame(rbind(imprinted_clusters['Cluster1', ],
                                          colMeans(imprinted_clusters[c('Cluster2a', 'Cluster2b', 'Cluster2c'), ])))
row.names(imprinted_clusters) <- c('Naive', 'Memory')
joint_names <- intersect(row.names(imprinted_bulk), colnames(imprinted_clusters))
imprinted_bulk <- imprinted_bulk[joint_names, ]
imprinted_clusters <- imprinted_clusters[, joint_names]
to_add <- data.frame(naiveB=imprinted_bulk$ImprintedBulk, memoryB=imprinted_bulk$ImprintedBulk, t(imprinted_clusters), Amplicon=joint_names)
to_add <- reshape2::melt(to_add, id=c('Amplicon', 'naiveB', 'memoryB'))
colnames(to_add)[4:5] <- c('CellType', 'Methylation')
to_add <- reshape2::melt(to_add, id=c('Amplicon', 'CellType', 'Methylation'))
colnames(to_add)[4:5] <- c('Bulk', 'BulkMethylation')
to_add$Correlation <- NA
to_add$CorrelationPVal <- NA
to_add$Type <- 'Imprinted'
to_plot <- as.data.frame(rbind(to_plot, to_add))
color_map <- c('1'='#fcbd7e',
               '2'='#fc3262',
               '3'='gray20',
               '4'=NA)
plot <- ggplot(to_plot, aes_string(x='Methylation', y='BulkMethylation'))+geom_point(aes(color=ifelse(Type=='Imprinted', '3', CellType)))+geom_smooth(aes(color=ifelse(Type=='Imprinted', '4', CellType)), method='lm',se=FALSE)+xlim(0, 1)+ylim(0,1)+
  facet_wrap(Bulk~CellType, nrow=2)+
  geom_text(aes(label=paste('Pearsons R²:', format(Correlation, digits=2))), x=0.35, y=0.9, check_overlap=TRUE, color='black')+
  plot_theme+
  scale_color_manual(values=color_map)+xlab('Pseudo-bulk methylation')+ylab('Bulk methylation')
ggsave(file.path(plot_path, 'F1_E_bulk_vs_singlecell.pdf'), plot, width=110, height=120, units='mm')
