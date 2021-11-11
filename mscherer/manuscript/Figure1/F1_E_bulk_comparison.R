################### F1_E_bulk_comparison.R ################### 
#'This file compares the pseudo-bulk obtained from averaging the methylation calls
#'in the clusters defined in the heatmap with the bulk values measured in
#'doi.org/10.1038/ng.3291


library(ggplot2)
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=15),
                    axis.text=element_text(color='black',size=15),
                    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    legend.key=element_rect(color='black', fill=NA))
color_map <- c('naive'='#fcbd7e',
               'memory1'='#fc6571',
               'memory2'='#fc3262')

filtered.counts <- read.table("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/tsv/BCells_Sample7_70_percent_good_performance.barcode.cell.distribution_with_MCL.tsv", row.names = 1, header=T)
cell_metadata <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/cell_metadata.csv',
                          row.names=1)
amplicon.info <- read.csv("/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/selected_amplicons.csv", row.names = 1)

selected <- ifelse(filtered.counts[row.names(cell_metadata), row.names(amplicon.info)]>0, 1, 0)

means_clusters <- aggregate(selected, by=list(cell_metadata$CellType), mean)
bulk_table <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
means_bulk <- bulk_table[, c('NBC.mean', 'MBC.mean')]
colnames(means_bulk) <- c('naiveB', 'memoryB')
row.names(means_bulk) <- bulk_table$amplicon
means_bulk <- as.data.frame(t(means_bulk))
means_bulk$'Group.1' <- row.names(means_bulk)
to_plot <- rbind(means_clusters, means_bulk[, colnames(means_clusters)[colnames(means_clusters)%in%colnames(means_bulk)]])
row.names(to_plot) <- to_plot$Group.1
to_plot <- as.data.frame(t(to_plot[, -1]))
to_plot$Amplicon <- row.names(to_plot)
to_plot <- reshape2::melt(to_plot, id=c('Amplicon', 'naiveB', 'memoryB'))
colnames(to_plot)[4:5] <- c('CellType', 'Methylation')
to_plot <- reshape2::melt(to_plot, id=c('Amplicon', 'CellType', 'Methylation'))
colnames(to_plot)[4:5] <- c('Bulk', 'BulkMethylation')
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
cors <- t(as.data.frame(apply(to_plot, 1, function(x){
  ct_sc <- x['CellType']
  ct_bulk <- x['Bulk']
  meth <- subset(to_plot, CellType==ct_sc & Bulk==ct_bulk, select=c('Methylation', 'BulkMethylation'))
  c(cor(meth$Methylation, meth$BulkMethylation),
    cor.test(meth$Methylation, meth$BulkMethylation)$p.value)
})))
to_plot$Correlation <- cors[,1]
to_plot$CorrelationPVal <- cors[,2]
plot <- ggplot(to_plot, aes_string(x='BulkMethylation', y='Methylation', color='CellType'))+geom_point()+geom_smooth(method='lm',se=FALSE)+xlim(0, 1)+ylim(0,1)+
  facet_grid(Bulk~CellType)+
  geom_text(aes(label=paste('Pearsons R²:', format(Correlation, digits=2))), x=0.5, y=0.9, check_overlap=TRUE, color='black')+
  plot_theme+
  scale_color_manual(values=color_map)+ylab('Bulk methylation')+xlab('Pseudo-bulk methylation')
ggsave(file.path(plot_path, 'F1_E_bulk_vs_singlecell.pdf'), plot, width=200, height=150, units='mm')
