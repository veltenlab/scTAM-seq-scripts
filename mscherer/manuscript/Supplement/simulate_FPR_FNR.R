set.seed(42)
compute_statistics <- function(observed_labels,
                               observed_bulk,
                               simulated_data){
  require(fpc)
  pseudo_bulks <- aggregate(simulated_data, by=list(observed_labels$CellType), mean)
  row.names(pseudo_bulks) <- pseudo_bulks[, 1]
  pseudo_bulks <- pseudo_bulks[, -1]
  pseudo_bulks <- as.data.frame(t(pseudo_bulks))[row.names(observed_bulk), colnames(observed_bulk)]
  avg_cor <- mean(diag(cor(pseudo_bulks, observed_bulk)))
  clust <- cutree(hclust(dist(simulated_data, 'binary'),
                         'ward.D2'),
                  3)
  j_1 <- sum(clust==c('HSC'=1,
                      'naive'=2,
                      'memory'=3)[observed_labels$CellType])/nrow(observed_labels)
  j_2 <- sum(clust==c('HSC'=1,
                      'naive'=3,
                      'memory'=2)[observed_labels$CellType])/nrow(observed_labels)
  j_3 <- sum(clust==c('HSC'=2,
                      'naive'=1,
                      'memory'=3)[observed_labels$CellType])/nrow(observed_labels)
  j_4 <- sum(clust==c('HSC'=2,
                      'naive'=3,
                      'memory'=1)[observed_labels$CellType])/nrow(observed_labels)
  j_5 <- sum(clust==c('HSC'=3,
                      'naive'=1,
                      'memory'=2)[observed_labels$CellType])/nrow(observed_labels)
  j_6 <- sum(clust==c('HSC'=3,
                      'naive'=2,
                      'memory'=1)[observed_labels$CellType])/nrow(observed_labels)
  jaccard <- max(c(j_1, j_2, j_3, j_4, j_5, j_6))
  return(c(Cor=avg_cor,
         Hamming=jaccard))
}

ampli_info <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/missionbio/BM/Sample11_70_percent_good_performance/tsv/colinfo.csv',
                       row.names = 1)
sel_data <- ampli_info[sample(1:nrow(ampli_info), 100), c('S1.mean',
                                                          'NBC.mean',
                                                          'MBC.mean')]
colnames(sel_data) <- c('HSC', 'naive', 'memory')
simulated_data <- apply(sel_data, c(1,2), function(x){
  rbinom(500, 1, x)
})
simulated_data <- rbind(simulated_data[,,1], rbind(simulated_data[,,2], simulated_data[,,3]))
cell_names <-  c(paste0('HSC', 1:500),
                 paste0('naive', 1:500),
                 paste0('memory', 1:500))
row.names(simulated_data) <- cell_names
annotation <- data.frame(CellType=c(rep('HSC', 500),
                                    rep('naive', 500),
                                    rep('memory', 500)))
row.names(annotation) <- cell_names
fprs <- fnrs <- seq(0, 0.5, by=0.025)
res <- c()
for(fpr in fprs){
  for(fnr in fnrs){
    print(paste('FPR:', fpr, 'FNR:', fnr))
    data_with_noise <- simulated_data
    if(fpr>0){
      sel_entries <- which(simulated_data==0, arr.ind = TRUE)
      sel_entries <- sel_entries[sample(1:nrow(sel_entries), nrow(sel_entries)*(fpr)), ]
      for(i in 1:nrow(sel_entries)){
        data_with_noise[sel_entries[i, 'row'], sel_entries[i, 'col']] <- 1
      }
    }
    if(fnr>0){
      sel_entries <- which(simulated_data>0, arr.ind = TRUE)
      sel_entries <- sel_entries[sample(1:nrow(sel_entries), nrow(sel_entries)*(fnr)), ]
      for(i in 1:nrow(sel_entries)){
        data_with_noise[sel_entries[i, 'row'], sel_entries[i, 'col']] <- 0
      }
    }
    stats <- compute_statistics(annotation,
                       sel_data,
                       data_with_noise)
    res <- rbind(res, c(fpr, fnr, stats))
  }
}
colnames(res) <- c('FPR', 'FNR', 'Correlation', 'HammingSimilarity')
to_plot <- reshape2::melt(as.data.frame(res), id=c('FPR', 'FNR'))
library(ggplot2)
library(ggsci)
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=10),
                    axis.text=element_text(color='black',size=8),
                    axis.ticks=element_line(color='black', size=.1),
                    strip.background = element_blank(),
                    legend.key=element_rect(color='black', fill=NA),
                    legend.title = element_blank(),
                    legend.text = element_text(color='black',size=8),
                    plot.title=element_blank())
to_plot$FNR <- as.factor(to_plot$FNR)
plot <- ggplot(to_plot, aes(x=FPR, y=value, color=FNR, fill=FNR))+geom_point()+geom_line()+facet_grid(variable~FNR)+plot_theme

to_plot$FNR <- as.numeric(as.character(to_plot$FNR))
to_plot$Product <- to_plot$FPR*to_plot$FNR
to_plot$Special <- 'FNR<=0.25'
to_plot[to_plot$FNR==0.075&to_plot$FPR==0.025, 'Special'] <- 'scTAM-seq'
to_plot[to_plot$FNR>0.25, 'Special'] <- 'FNR>0.25'
to_plot$Special <- factor(to_plot$Special,
                          levels=c('FNR>0.25', 'FNR<=0.25', 'scTAM-seq'))
to_plot$FNR <- as.factor(to_plot$FNR)
to_plot$FPR <- as.factor(to_plot$FPR)
to_plot <- to_plot[order(to_plot$Special), ]
plot <- ggplot(to_plot, aes(x=log2(Product), y=value, color=Special, size=Special))+geom_smooth(se=FALSE)+geom_point()+facet_grid(variable~.)+plot_theme+
  xlab('log2(FPRxFNR)')+ylab("")+scale_color_jama()+scale_size_manual(values=c('FNR>0.25'=1,
                                                                        'FNR<=0.25'=1,
                                                                        'scTAM-seq'=3),
                                                                      guide='none')
ggsave('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/revision/FPR_FNR_simulation.pdf',
    width = 150, height = 100, unit='mm')
