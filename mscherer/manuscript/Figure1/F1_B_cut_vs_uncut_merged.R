############## F1_B_cut_vs_uncut.R ##############
#' This file generates the plot comparing the dropout rate per amplicon (type) in the cut versus the uncut
#' sample. Dropout is defined as those cells having 0 reads at this particular amplicon

library(ggplot2)
cut <- 'Sample8_70_percent_good_performance'
s7 <- 'Sample7_70_percent_good_performance'
uncut <- 'Sample6_70_percent_good_performance'
plot.path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/Sample8/'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=6),
               axis.text=element_text(color='black',size=5),
               axis.ticks=element_line(color='black', size=.25),
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
dat_cut <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/", cut, "/tsv/", cut, ".barcode.cell.distribution.tsv"), 
                  sep="\t", 
                  header=T)
#doublet <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/old/Figure1/reclustering_doublet_names.csv')
#double2 <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/tsv/doublet_scores_DoubletDetection.csv')
#doublets <- c(doublet$x, double2$Barcode[double2$DoubletDetectionLabel==1])
#doublets <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/' ,
#  cut, '/tsv/doublet_scores_DoubletDetection.csv'))
#doublets <- doublets$Barcode[doublets$DoubletDetectionLabel==1]
rowinfo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/' ,
                           cut, '/tsv/rowinfo.csv'),
                    row.names = 1)
dat_cut <- dat_cut[row.names(rowinfo), ]
dat_s7 <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/", s7, "/tsv/", s7, ".barcode.cell.distribution.tsv"), 
                      sep="\t", 
                      header=T)
doublets <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/' ,
                            s7, '/tsv/doublet_scores_DoubletDetection.csv'))
doublets <- doublets$Barcode[doublets$DoubletDetectionLabel==1]
dat_s7 <- dat_s7[!(row.names(dat_s7)%in%doublets), ]
colliding.barcodes <- plyr::count(c(row.names(dat_cut), row.names(dat_s7)))
colliding.barcodes <- colliding.barcodes[colliding.barcodes$freq>1, 'x']
dat_cut <- dat_cut[-(which(row.names(dat_cut)%in%colliding.barcodes)), ]
dat_s7 <- dat_s7[-(which(row.names(dat_s7)%in%colliding.barcodes)), ]
#dat_cut <- as.data.frame(rbind(dat_cut, dat_s7))
dat_uncut <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/",
                               uncut, "/tsv/", uncut, ".barcode.cell.distribution.tsv"), 
                      sep="\t", 
                      header=T)
doublet <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', 
  uncut, '/tsv/doublet_scores_DoubletDetection.csv'))
doublets <- c(doublet$Barcode[doublet$DoubletDetectionLabel==1])
dat_uncut <- dat_uncut[!(row.names(dat_uncut)%in%doublets), ]
ampli_info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
in_all <- intersect(colnames(dat_cut), intersect(colnames(dat_uncut), row.names(ampli_info)))
dropout_cut <- apply(dat_cut[, in_all], 2, function(x){
  1-(sum(x==0)/length(x))
})
dropout_uncut <- apply(dat_uncut[, in_all], 2, function(x){
  1-(sum(x==0)/length(x))
})
ampli_info <- ampli_info[in_all, ]
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
to_plot <- data.frame(Cut=dropout_cut,
                      Uncut=dropout_uncut,
                      Type=type_clear[ampli_info$Type.of.amplicon])
plot_theme <- plot_theme+theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5))
plot <- ggplot(to_plot, aes(x=Uncut, y=Cut, color=Type))+
  geom_point()+geom_abline(slope=1, intercept=0)+facet_wrap(Type~.)+
  labs(x='Fraction of cells with reads in undigested sample', y='Fraction of cells with reads in digested sample', col='Amplicon Type')+  plot_theme+scale_color_manual(values=colors_amplicons)
ggsave(file.path(plot.path, 'F1_B_cut_vs_uncut_all.pdf'), plot, width=175, height=100, units='mm')
to_plot <- to_plot[to_plot$Type%in%c("No HhaI cutsite",
                                     "Differential CpG Bcells",
                                     "Always unmethylated CpG in Bcells",
                                     "Always methylated CpG in Bcells"),
                   ]
to_plot$Type <- factor(to_plot$Type, levels=c("No HhaI cutsite",
                                              "Always methylated CpG in Bcells",
                                              "Always unmethylated CpG in Bcells", 
                                              "Differential CpG Bcells"))
plot <- ggplot(to_plot, aes(x=Uncut, y=Cut, color=Type))+
  geom_point(size=.75)+geom_abline(slope=1, intercept=0, size=.25)+facet_wrap(Type~., nrow = 1)+
  labs(x='Fraction of cells with reads in undigested sample', y='Fraction of cells with reads in digested sample', col='Amplicon Type')+
plot_theme+scale_color_manual(values=colors_amplicons)
ggsave(file.path(plot.path, 'F1_B_cut_vs_uncut_selected.pdf'), plot, width=110, height=33, units='mm')
