library(ggplot2)
library(rstan)
plot_path <- '/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure1/Sample8/correction/'
cut <- 'Sample8_70_percent_good_performance'
uncut <- 'Sample6_70_percent_good_performance'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_line(color='black', size=.25),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')
color_map <- c('Naive'='#fcbd7e',
               'Memory'='#fc00bd')

filtered.counts <- read.table(paste0("/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/", cut, "/tsv/", cut, ".barcode.cell.distribution.tsv"), 
                              sep="\t", 
                              header=T)
cell_metadata <- read.csv(paste0('/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/', cut, '/tsv/rowinfo.csv'),
                          row.names = 1)
all.amplicons <- read.table('/users/lvelten/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
# Only using the best-performing amplicons
#sel.amplicons <- row.names(amplicon.info)
# Using all amplicons
sel.amplicons <- row.names(all.amplicons[grepl('CpG.B.cell.diff', all.amplicons$Type.of.amplicon),])
non_cut <- read.table(paste0("/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/", uncut, "/tsv/", uncut, ".barcode.cell.distribution.tsv"), 
                      sep="\t", 
                      header=T)
non_cut_doublets <- read.csv(paste0('/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/', 
                                    uncut, '/tsv/doublet_scores_DoubletDetection.csv'))
non_cut <- non_cut[non_cut_doublets[which(non_cut_doublets$DoubletDetectionLabel==0), 'Barcode'], ]
allelic_dropout <- apply(non_cut[, sel.amplicons], 2, function(x){
  sum(x==0)/length(x)
})
sel.amplicons <- sel.amplicons[allelic_dropout>0.1]
allelic_dropout <- allelic_dropout[sel.amplicons]

selected <- ifelse(filtered.counts[row.names(cell_metadata), sel.amplicons]>0, 1, 0)

cell_metadata$CellType <- ifelse(cell_metadata$CellType_broad%in%'naive B-cells', 'naive B-cells', 'memory B-cells')
vals_naive <- apply(selected[cell_metadata$CellType_broad%in%'naive B-cells', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})
vals_memory <- apply(selected[cell_metadata$CellType_broad%in%'memory B-cells', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})
corrected_values <- sapply(sel.amplicons, function(ampli){
  print(ampli)
  data_naive <- list(n0=vals_naive[1, ampli], n1=vals_naive[2, ampli], p=allelic_dropout[ampli])
  sm_naive <- stan("basic.stan", data = data_naive)
  data_memory <- list(n0=vals_memory[1, ampli], n1=vals_memory[2, ampli], p=allelic_dropout[ampli])
  sm_memory <- stan("basic.stan", data = data_memory)
  c(Naive=get_posterior_mean(sm_naive)['m', 'mean-all chains'],
    Memory=get_posterior_mean(sm_memory)['m', 'mean-all chains'])
})

bulk_table <- read.table('/users/lvelten/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
means_bulk <- bulk_table[, c('NBC.mean', 'MBC.mean')]
colnames(means_bulk) <- c('naiveB', 'memoryB')
row.names(means_bulk) <- bulk_table$amplicon
means_bulk <- as.data.frame(t(means_bulk))
to_plot <- rbind(corrected_values, means_bulk[, colnames(corrected_values)[colnames(corrected_values)%in%colnames(means_bulk)]])
to_plot <- as.data.frame(t(to_plot))
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
#ggsave(file.path(plot_path, 'F1_E_bulk_vs_singlecell_corrected_stan.pdf'), plot, width=200, height=150, units='mm')
