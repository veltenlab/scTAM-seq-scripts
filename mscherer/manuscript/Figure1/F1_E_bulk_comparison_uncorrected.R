################### F1_E_bulk_comparison_uncorrected.R ################### 
#'This file compares the pseudo-bulk obtained by computing the number of methylated
#'cells per cluster with the bulk values from:
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
color_map <- c('naive B-cells'='#fcbd7e',
               'memory B-cells'='#fc3262')
sample <- 'Sample7_70_percent_good_performance'

dat_cut <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/", sample, "/tsv/", cut, ".barcode.cell.distribution.tsv"), 
                      sep="\t", 
                      header=T)
rowinfo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/rowinfo.csv'),
                    row.names = 1)
ampli_info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
sel_amplis <- row.names(ampli_info)[grepl('CpG.B.cell.diff', ampli_info$Type.of.amplicon)]
dat_cut <- dat_cut[row.names(rowinfo), sel_amplis]
means_clusters <- aggregate(dat_cut, by=list(rowinfo$CellType_broad), function(x){
  sum(x>0)/length(x)
})
c.names <- means_clusters[, 1]
means_clusters <- as.data.frame(t(means_clusters[, -1]))
colnames(means_clusters) <- c.names
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
plot <- ggplot(to_plot, aes_string(x='BulkMethylation', y='Methylation', color='CellType'))+geom_point()+geom_smooth(method='lm',se=FALSE)+xlim(0, 1)+ylim(0,1)+
  facet_wrap(Bulk~CellType, nrow=2)+
  geom_text(aes(label=paste('Pearsons R²:', format(Correlation, digits=2))), x=0.35, y=0.9, check_overlap=TRUE, color='black')+
  plot_theme+
  scale_color_manual(values=color_map)+ylab('Pseudo-bulk methylation')+xlab('Bulk methylation')
ggsave(file.path(plot_path, 'F1_E_bulk_vs_singlecell_uncorrected.pdf'), plot, width=110, height=120, units='mm')
