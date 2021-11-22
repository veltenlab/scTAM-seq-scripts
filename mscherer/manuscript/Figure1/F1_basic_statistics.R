############## F1_basic_statistics.R ##############
#' This file generates basic information about the experiment, including false positive rate, false
#' negative rate and dropout rates

library(ggplot2)
cut <- 'BCells_Sample7_70_percent_good_performance'
uncut <- 'BCells_Sample6_70_percent_good_performance'
plot.path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=10),
                    axis.text=element_text(color='black',size=8),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')
type_clear <- c("mut.only"="Mutation only",
                "NonHhaI"="No HhaI cutsite",
                "CpG.B.cell.diff"="Differential CpG Bcells",
                "CpG.MCL"="MCL-specific CpG",
                "CpG.MCL.multiple.cutsites"="MCL-specific CpG",
                "CpG.always.unmeth.B"="Always unmethylated CpG in Bcells",
                "CpG.B.cell.diff.and.MCL"="Differential CpG Bcells AND MCL-specific CpG",        
                "CpG.always.meth.B"="Always methylated CpG in Bcells",
                "CpG.imprinted.multiple.cutsites"="Imprinted CpG multiple",
                "CpG.imprinted"="Imprinted CpG")
colors_amplicons <- c("Mutation only"="#e5c494",
                      "Differential CpG Bcells"="#fc8d62",
                      "MCL-specific CpG"="#ffd92f",
                      "Differential CpG Bcells AND MCL-specific CpG"="#fc8d62",
                      "Always unmethylated CpG in Bcells"="#a6d854",
                      "Always methylated CpG in Bcells"="#66c2a5",
                      "Imprinted CpG"="#e78ac3",
                      "Imprinted CpG multiple"="#e78ac3",
                      "No HhaI cutsite"="#b3b3b3")
dat_cut <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/", cut, "/tsv/", cut, ".barcode.cell.distribution_with_MCL.tsv"), 
                      sep="\t", 
                      header=T)
rowinfo <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/rowinfo_reclustering_doublet.csv',
                    row.names = 1)
doublet <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/reclustering_doublet_names.csv')
double2 <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/tsv/doublet_scores_DoubletDetection.csv')
doublets <- c(doublet$x, double2$Barcode[double2$DoubletDetectionLabel==1])
dat_cut <- dat_cut[!(row.names(dat_cut)%in%doublets), ]
dat_uncut <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/", uncut, "/tsv/", uncut, ".barcode.cell.distribution.tsv"), 
                        sep="\t", 
                        header=T)
doublet <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample6_70_percent_good_performance/tsv/doublet_scores_DoubletDetection.csv')
doublets <- c(doublet$Barcode[doublet$DoubletDetectionLabel==1])
dat_uncut <- dat_uncut[!(row.names(dat_uncut)%in%doublets), ]
ampli_info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
selected_amplicons <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/selected_amplicons.csv',
                               row.names = 1)
in_all <- intersect(colnames(dat_cut), intersect(colnames(dat_uncut), row.names(ampli_info)))
ampli_info <- ampli_info[][in_all, ]

fpr <- ifelse(dat_cut[, ampli_info$Type.of.amplicon%in%"CpG.always.unmeth.B"]>0, 1, 0)
fpr <- colMeans(fpr)
fnr <- ifelse(dat_cut[, ampli_info$Type.of.amplicon%in%"CpG.always.meth.B"]>0, 0, 1)
fnr <- colMeans(fnr)
dropout <- ifelse(dat_uncut[, row.names(selected_amplicons)]>0, 0, 1)
dropout <- colMeans(dropout)
to_plot <- data.frame(Type=c(rep('FPR', length(fpr)),
                             rep('FNR', length(fnr)),
                             rep('Dropout', length(dropout))),
                     Value=c(fpr, fnr, dropout))
plot <- ggplot(to_plot, aes(x=Type, y=Value))+geom_boxplot(color='black', fill='gray80')+plot_theme
ggsave('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/F1_basic_statistics.pdf', 
       plot,
       width=50,
       height=35,
       units='mm')
