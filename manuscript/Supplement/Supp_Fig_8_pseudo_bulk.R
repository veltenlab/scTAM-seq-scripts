################### Supp_Fig_6.R ################### 
#'This file compares the pseudo-bulk obtained from averaging the methylation calls
#'in the clusters defined in the heatmap with the bulk values measured in
#'doi.org/10.1038/ng.3291


library(ggplot2)
plot_path <- '~'
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
sample <- 'GSM5935921_BM_HhaI'
means_clusters <- read.csv(paste0('../../dropout_modeling/', sample, '/corrected_values.csv'),
                           row.names=1)
colinfo <- read.table("../../misc/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt", header=T,
                      row.names=8)
colinfo <- colinfo[, !(colnames(colinfo)%in%c('GC.mean', 'PB.mean'))]
colnames(colinfo) <- c('S1.mean'='Progenitors',
                       'S2.mean'='PreB',
                       'S3.mean'='Pre-BII',
                       'S4.mean'='ImmatureB',
                       'NBC.mean'='naiveB',
                       'MBC.mean'='memoryB')[colnames(colinfo)]
colinfo$ImmatureB <- rowMeans(colinfo[, c('Pre-BII', 'ImmatureB')])
to_plot <- data.frame(means_clusters, colinfo[row.names(means_clusters), c('Progenitors', 'PreB', 'ImmatureB', 'naiveB', 'memoryB')])
to_plot$Amplicon <- row.names(to_plot)
to_plot <- reshape2::melt(to_plot, id=c('Amplicon', 'Progenitors', 'PreB', 'ImmatureB', 'naiveB', 'memoryB'))
colnames(to_plot)[7:8] <- c('SingleCell', 'Methylation')
to_plot <- reshape2::melt(to_plot, id=c('Amplicon', 'SingleCell', 'Methylation'))
colnames(to_plot)[4:5] <- c('Bulk', 'BulkMethylation')
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
color_map <- c('S1.cells'='#bdfc7e',
               'S2.cells'='#d7ef7e',
               'S3.S4.cells'='#fcef7e',
               'naive.B.cells'='#fcbd7e',
               'memory.B.cells'='#fc6571')
to_plot$SingleCell <- factor(to_plot$SingleCell, levels=c('S1.cells', 'S2.cells', 'S3.S4.cells', 'naive.B.cells', 'memory.B.cells'))
to_plot$Bulk <- factor(to_plot$Bulk, levels=c('Progenitors', 'PreB', 'ImmatureB', 'naiveB', 'memoryB'))
plot <- ggplot(to_plot, aes_string(x='BulkMethylation', y='Methylation', color='SingleCell'))+geom_point(size=.1)+geom_smooth(method='lm',se=FALSE, size=.5)+xlim(0, 1)+ylim(0,1)+
  facet_wrap(SingleCell~Bulk, nrow=5)+
  geom_text(aes(label=paste('R²:', format(Correlation, digits=2))), x=0.2, y=0.9, check_overlap=TRUE, color='black', size=1.75, fontface = "bold")+
  plot_theme+
  scale_color_manual(values=color_map)+ylab('Pseudo-bulk methylation')+xlab('Bulk methylation')
ggsave(file.path(plot_path, 'Supp_Fig_6.pdf'), plot, width=100, height=100, units='mm')

