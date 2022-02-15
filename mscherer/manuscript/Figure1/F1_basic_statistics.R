############## F1_basic_statistics.R ##############
#' This file generates basic information about the experiment, including false positive rate, false
#' negative rate and dropout rates

library(ggplot2)
cut <- 'Sample11_70_percent_good_performance'
uncut <- 'Sample12_70_percent_good_performance'
plot.path <- '/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure1/Sample11/'
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
dat_cut <- read.table(paste0("/users/lvelten/project/Methylome/analysis/missionbio/BM/", cut, "/tsv/", cut, ".barcode.cell.distribution.tsv"), 
                      sep="\t", 
                      header=T)
rowinfo <- read.csv(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/', cut, '/tsv/rowinfo.csv'),
                    row.names = 1)
dat_uncut <- read.table(paste0("/users/lvelten/project/Methylome/analysis/missionbio/BM/", uncut, "/tsv/", uncut, ".barcode.cell.distribution.tsv"), 
                        sep="\t", 
                        header=T)
rowinfo_uncut <- read.csv(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/', uncut, '/tsv/doublet_scores_DoubletDetection.csv'),
                    row.names = 2)

doublet <- read.csv(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/', 
                           uncut, '/tsv/doublet_scores_DoubletDetection.csv'))
doublets <- c(doublet$Barcode[which(doublet$DoubletDetectionLabel==1)])
dat_uncut <- dat_uncut[!(row.names(dat_uncut)%in%doublets), ]
ampli_info <- read.table('/users/lvelten/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
selected_amplicons <- read.table(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/', 
                                      uncut, '/tsv/selected_amplicons.tsv'),
                               row.names = 1)
in_all <- intersect(colnames(dat_cut), intersect(colnames(dat_uncut), row.names(ampli_info)))
ampli_info <- ampli_info[][in_all, ]

fpr <- apply(dat_cut[row.names(rowinfo), row.names(ampli_info)[grepl("CpG.always.unmeth.B", ampli_info$Type.of.amplicon)]], 2, function(x){
  sum(x>0)/length(x)
})
fnr <- apply(dat_uncut[row.names(rowinfo_uncut[which(rowinfo_uncut$DoubletDetectionLabel==0), ]), row.names(ampli_info)[grepl("CpG.B.cell.diff", ampli_info$Type.of.amplicon)]], 2, function(x){
  sum(x==0)/length(x)
})
#fpr <- ifelse(dat_cut[row.names(rowinfo), ampli_info$Type.of.amplicon%in%"CpG.always.unmeth.B"]>0, 1, 0)
#fpr <- colMeans(fpr)
#fnr <- ifelse(dat_uncut[row.names(rowinfo_uncut[which(rowinfo_uncut$DoubletDetectionLabel==0), ]), row.names(ampli_info)[ampli_info$Type.of.amplicon%in%"CpG.B.cell.diff"]]>0, 0, 1)
# fnr <- ifelse(dat_uncut[row.names(rowinfo_uncut[which(rowinfo_uncut$DoubletDetectionLabel==0), ]), ampli_info$Type.of.amplicon%in%c("CpG.B.cell.diff",
#                                                                                                                                     "NonHhaI",
#                                                                                                                                     "CpG.always.unmeth.B",
#                                                                                                                                     "CpG.always.meth.B")]>0, 0, 1)
#fnr <- colMeans(fnr)
#dropout <- ifelse(dat_cut[row.names(rowinfo), ampli_info$Type.of.amplicon%in%"NonHhaI"]>0, 0, 1)
#dropout <- colMeans(dropout)
to_plot <- data.frame(Type=c(rep('FPR', length(fpr)),
                             rep('FNR', length(fnr))),
                     Value=c(fpr, fnr))
plot <- ggplot(to_plot, aes(x=Type, y=Value))+geom_boxplot(color='black', fill='gray80', size=.25, outlier.size=.25)+plot_theme+
  xlab('')+ylab('')
ggsave('/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure1/Sample11/F1_basic_statistics.pdf', 
       plot,
       width=40,
       height=40,
       units='mm')
