library(ggplot2)
plot_path <- '~'
sample <- 'GSM5935921_BM_HhaI'
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
means_clusters <- read.csv(file.path('../../dropout_modeling', sample, 'corrected_values.csv'),
                           row.names=1)
colinfo <- read.table('../../misc/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt',
                      header=T,
                      row.names=8)
special <- read.csv(file.path('../../misc/', sample, 
                              'tsv/colinfo_BM.csv'),
                    row.names=1)
colinfo[row.names(special), 'Special'] <- special$Special
colinfo$Special[is.na(colinfo$Special)] <- 'Normal'
colinfo <- colinfo[, !(colnames(colinfo)%in%c('GC.mean', 'PB.mean'))]
colnames(colinfo) <- c('S1.mean'='Progenitors',
                       'S2.mean'='PreB',
                       'S3.mean'='Pre-BII',
                       'S4.mean'='ImmatureB',
                       'NBC.mean'='naiveB',
                       'MBC.mean'='memoryB',
                       'Special'='Special')[colnames(colinfo)]
colinfo$ImmatureB <- rowMeans(colinfo[, c('Pre-BII', 'ImmatureB')])
to_plot <- data.frame(means_clusters, colinfo[row.names(means_clusters), c('Progenitors', 'PreB', 'ImmatureB', 'naiveB', 'memoryB', 'Special')])
to_plot$Amplicon <- row.names(to_plot)
to_plot <- reshape2::melt(to_plot, id=c('Amplicon', 'Progenitors', 'PreB', 'ImmatureB', 'naiveB', 'memoryB', 'Special'))
colnames(to_plot)[8:9] <- c('SingleCell', 'Methylation')
to_plot <- reshape2::melt(to_plot, id=c('Amplicon', 'SingleCell', 'Methylation', 'Special'))
colnames(to_plot)[5:6] <- c('Bulk', 'BulkMethylation')
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
to_plot$BulkMethylation <- as.numeric(to_plot$BulkMethylation)
to_plot <- na.omit(to_plot)
cors <- t(as.data.frame(apply(to_plot, 1, function(x){
  ct_sc <- x['SingleCell']
  ct_bulk <- x['Bulk']
  meth <- subset(to_plot, SingleCell==ct_sc & Bulk==ct_bulk, select=c('Methylation', 'BulkMethylation'))
  c(cor(meth$Methylation, meth$BulkMethylation),
    cor.test(meth$Methylation, meth$BulkMethylation)$p.value)
})))
to_plot$Correlation <- cors[,1]
to_plot$CorrelationPVal <- cors[,2]
to_plot$Bulk <- factor(to_plot$Bulk, levels=c('Progenitors', 'PreB', 'ImmatureB', 'naiveB', 'memoryB'))
color_map <- c('S1cells'='#bdfc7e',
               'S2cells'='#d7ef7e',
               'S3S4cells'='#fcef7e',
               'NaiveB'='#fcbd7e',
               'MemoryB'='#fc6571')
to_plot$SingleCell <- factor(to_plot$SingleCell, levels=c('S1cells', 'S2cells', 'S3S4cells', 'NaiveB', 'MemoryB'))
to_plot$Bulk <- factor(to_plot$Bulk, levels=c('Progenitors', 'PreB', 'ImmatureB', 'naiveB', 'memoryB'))
plot <- ggplot(to_plot, aes_string(x='BulkMethylation', y='Methylation'))+geom_point(aes(color=ifelse(Special=='Normal', 'gray10', 'black')), size=.1)+geom_smooth(method='lm',se=FALSE, size=.5, color='#808080')+xlim(0, 1)+ylim(0,1)+
  facet_wrap(SingleCell~Bulk, nrow=5)+
  geom_text(aes(label=paste('r:', format(Correlation, digits=2))), x=0.2, y=0.9, check_overlap=TRUE, color='black', size=1.75, fontface = "bold")+
  plot_theme+ylab('Pseudo-bulk methylation')+xlab('Bulk methylation')
ggsave(file.path(plot_path, 'F2_A_bulk_vs_singlecell.pdf'), plot, width=100, height=100, units='mm')

