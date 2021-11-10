library(ggplot2)
treated_sample <- read.table('/users/lvelten/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/tsv/BCells_Sample7_70_percent_good_performance.barcode.cell.distribution_with_MCL.tsv')
untreated_sample <- read.table('/users/lvelten/project/Methylome/analysis/missionbio/tapestri/BCells_Sample6_70_percent_good_performance/tsv/BCells_Sample6_70_percent_good_performance.barcode.cell.distribution.tsv')
ampli.info <- read.table('/users/lvelten/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
sel.amplis <- row.names(ampli.info)[ampli.info$'Type.of.amplicon'%in%c('CpG.always.meth.B', 'NonHhaI')]
dropout_amplicons_untreated <- apply(untreated_sample[, sel.amplis], 2, function(x){
  sum(x==0)/length(x)
})
dropout_amplicons_treated <- apply(treated_sample[, sel.amplis], 2, function(x){
  sum(x==0)/length(x)
})
sel.amplis <- sel.amplis[dropout_amplicons_untreated<=0.8&
                         dropout_amplicons_untreated>=0.15&
                         dropout_amplicons_treated<=0.8]#&
#                         dropout_amplicons_treated>=0.15]
dropout_treated <- apply(treated_sample[, sel.amplis], 1, function(x){
  sum(x==0)/length(x)
})
dropout_untreated <- apply(untreated_sample[, sel.amplis], 1, function(x){
  sum(x==0)/length(x)
})
#dropout_treated <- dropout_treated[dropout_treated>0.05]
to_plot <- data.frame(Sample=factor(c(rep('Treated', length(dropout_treated)),
                               rep('Untreated', length(dropout_untreated))), levels=c('Untreated', 'Treated')),
                      Dropout=c(dropout_treated, dropout_untreated))
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
plot <- ggplot(to_plot, aes(x=Dropout, y=..count.., color=Sample))+geom_histogram(fill=NA, position='identity', binwidth=0.005)+theme

library(diptest)
dip.test(dropout_treated)
library(multimode)
modetest(dropout_treated)
res <- locmodes(dropout_treated, mod0=2, display=TRUE)
mode <- res$locations[1]
do.ass <- ifelse(dropout_treated<mode, 'Doublet', 'Singlet')
plyr::count(do.ass)
#do.ass <- ifelse(dropout_treated<0.25, 'Doublet', 'Singlet')
write.csv(do.ass, '/users/lvelten/project/Methylome/analysis/doublet_detection/Hha_dropout/doublets.csv')
