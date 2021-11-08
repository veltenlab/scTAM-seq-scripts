###################### plot_amplicons_in_undigested.R ###############################
#' This file creates boxplots for the reads for a number of specifies amplicons.

library(ggplot2)
library(reshape2)
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15),
               axis.text.x=element_text(angle=45, hjust=1))
sample <- 'BCells_Sample7_70_percent_good_performance'
amplicons <- c('AMPL131008',
               'AMPL131358',
               'AMPL130912',
               'AMPL131007',
               'AMPL131138',
               'AMPL131332',
               'AMPL131166',
               'AMPL131130',
               'AMPL131131',
               'AMPL131139',
               'AMPL131140')

input <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                    sep='\t',
                    header=TRUE)
to.plot <- input[, amplicons]
to.plot <- reshape2::melt(to.plot)
colnames(to.plot) <- c('Amplicon', 'ReadCount')
type <- rep('Intermediately Methylated',nrow(to.plot))
type[to.plot$Amplicon%in%c('AMPL131008', 'AMPL131358')] <- 'Fully Methylated'
type[to.plot$Amplicon%in%c('AMPL131130', 'AMPL131131', 'AMPL131139', 'AMPL131140')] <- 'Dropout'
to.plot$Type <- type
plot <- ggplot(to.plot,aes(x=Amplicon, y=ReadCount, fill=Type))+geom_violin()+geom_boxplot(alpha=0.25)+theme
