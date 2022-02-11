############## F1_B_cut_vs_uncut.R ##############
#' This file generates the plot comparing the dropout rate per amplicon (type) in the cut versus the uncut
#' sample. Dropout is defined as those cells having 0 reads at this particular amplicon

library(ggplot2)
library(gridExtra)
all_cut <- c("/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample7_70_percent_good_performance/tsv/Sample7_70_percent_good_performance.barcode.cell.distribution.tsv",
             "/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/tsv/Sample8_70_percent_good_performance.barcode.cell.distribution.tsv",
             "/users/lvelten/project/Methylome/analysis/missionbio/BM/Sample11_70_percent_good_performance/tsv/Sample11_70_percent_good_performance.barcode.cell.distribution.tsv")
all_uncut <- c("/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample6_70_percent_good_performance/tsv/Sample6_70_percent_good_performance.barcode.cell.distribution.tsv",
           "/users/lvelten/project/Methylome/analysis/missionbio/BM/Sample12_70_percent_good_performance/tsv/Sample12_70_percent_good_performance.barcode.cell.distribution.tsv")
cut_names <- c('Blood1', 'Blood2', 'BM1')
uncut_names <- c('BloodControl', 'BMControl')
cut_doublets <-  c("/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample7_70_percent_good_performance/tsv/all_doublets.csv",
                   "/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/tsv/doublet_scores_DoubletDetection.csv",
                   "/users/lvelten/project/Methylome/analysis/missionbio/BM/Sample11_70_percent_good_performance/tsv/doublet_scores_DoubletDetection.csv")
uncut_doublet <-  c("/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/Sample6_70_percent_good_performance/tsv/doublet_scores_DoubletDetection.csv",
                    "/users/lvelten/project/Methylome/analysis/missionbio/BM/Sample12_70_percent_good_performance/tsv/doublet_scores_DoubletDetection.csv")
plot.path <- '/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Supplement/'
ampli_info <- read.table('/users/lvelten/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
sel_amplis <- row.names(ampli_info[ampli_info$Type.of.amplicon%in%'CpG.always.meth.B', ])
dat_uncut_s6 <- read.table(all_uncut[1], 
                        sep="\t", 
                        header=T)
doublets_uncut_s6 <- read.csv(uncut_doublet[1])
doublets_uncut_s6 <- c(doublets_uncut_s6$Barcode[doublets_uncut_s6$DoubletDetectionLabel==1])
dat_uncut_s6 <- dat_uncut_s6[!(row.names(dat_uncut_s6)%in%doublets_uncut_s6), ]
sum_reads_s6 <- sum(rowSums(dat_uncut_s6))
dat_uncut_s12 <- read.table(all_uncut[2], 
                           sep="\t", 
                           header=T)
doublets_uncut_s12 <- read.csv(uncut_doublet[2])
doublets_uncut_s12 <- c(doublets_uncut_s12$Barcode[doublets_uncut_s12$DoubletDetectionLabel==1])
dat_uncut_s12 <- dat_uncut_s12[!(row.names(dat_uncut_s12)%in%doublets_uncut_s12), ]
sum_reads_s12 <- sum(rowSums(dat_uncut_s12))
dat_uncut_s12_sub_75 <- apply(dat_uncut_s12, c(1, 2), function(x){
  round(x*sample(seq(0.5, 1, by=0.01), 1))
})
sum_reads_s12_sub_75 <- sum(rowSums(dat_uncut_s12_sub_75))
dat_uncut_s12_sub_50 <- apply(dat_uncut_s12, c(1, 2), function(x){
  round(x*sample(seq(0, 1, by=0.01), 1))
})
sum_reads_s12_sub_50 <- sum(rowSums(dat_uncut_s12_sub_50))
dat_uncut_s12_sub_25 <- apply(dat_uncut_s12, c(1, 2), function(x){
  round(x*sample(seq(0, 0.5, by=0.01), 1))
})
sum_reads_s12_sub_25 <- sum(rowSums(dat_uncut_s12_sub_25))
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=6),
               axis.text=element_text(color='black',size=5),
               axis.ticks=element_line(color='black', size=.25),
               strip.background = element_blank(),
               strip.text.x = element_blank(),
               legend.key=element_rect(color=NA, fill=NA),
               legend.position='none')
dropout_uncut_s6 <- apply(dat_uncut_s6, 2, function(x){
  1-(sum(x==0)/length(x))
})
dropout_uncut_s12 <- apply(dat_uncut_s12, 2, function(x){
  1-(sum(x==0)/length(x))
})
dropout_uncut_s12_sub_75 <- apply(dat_uncut_s12_sub_75, 2, function(x){
  1-(sum(x==0)/length(x))
})
dropout_uncut_s12_sub_50 <- apply(dat_uncut_s12_sub_50, 2, function(x){
  1-(sum(x==0)/length(x))
})
dropout_uncut_s12_sub_25 <- apply(dat_uncut_s12_sub_25, 2, function(x){
  1-(sum(x==0)/length(x))
})
to_plot_dropout <- data.frame(Dropout=c(dropout_uncut_s6,
                                        dropout_uncut_s12,
                                        dropout_uncut_s12_sub_25,
                                        dropout_uncut_s12_sub_50,
                                        dropout_uncut_s12_sub_75),
                              Sample=factor(c(rep('Sample6',length(dropout_uncut_s6)),
                                       rep('Sample12',length(dropout_uncut_s6)),
                                       rep('Sample12_sub_25',length(dropout_uncut_s6)),
                                       rep('Sample12_sub_50',length(dropout_uncut_s6)),
                                       rep('Sample12_sub_75',length(dropout_uncut_s6))), levels=c('Sample12', 'Sample12_sub_75', 'Sample12_sub_50', 'Sample12_sub_25', 'Sample6')))
to_plot_reads <- data.frame(Reads=c(sum_reads_s6,
                                      sum_reads_s12,
                                      sum_reads_s12_sub_25,
                                      sum_reads_s12_sub_50,
                                      sum_reads_s12_sub_75),
                              Sample=factor(c('Sample6',
                                       'Sample12',
                                       'Sample12_sub_25',
                                       'Sample12_sub_50',
                                       'Sample12_sub_75'), levels=c('Sample12', 'Sample12_sub_75', 'Sample12_sub_50', 'Sample12_sub_25', 'Sample6')))
plot_dropout <- ggplot(to_plot_dropout, aes(x=Sample, y=1-Dropout))+geom_boxplot()+plot_theme+xlab('Dropout')
plot_reads <- ggplot(to_plot_reads, aes(x=Sample, y=Reads))+geom_histogram(stat='identity')+plot_theme+
  geom_text(aes(label=Reads), y=1e+8)
pdf(file.path(plot.path, 'dropout_reads_subsampling.pdf'), 
    width=6,
    height=6)
grid.arrange(plot_dropout, plot_reads)
dev.off()
