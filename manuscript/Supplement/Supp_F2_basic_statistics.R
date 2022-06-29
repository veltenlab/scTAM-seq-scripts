############## F1_basic_statistics.R ##############
#' This file generates basic information about the experiment, including false positive rate, false
#' negative rate and dropout rates

library(ggplot2)
cut <- 'GSM5935918_Blood_HhaI'
uncut <- 'GSM5935920_Blood_undigested'
plot.path <- '~'
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
dat_cut <- read.table(paste0("../../data/", cut, "/tsv/", cut, ".tsv.gz"), 
                      sep="\t", 
                      header=T)
rowinfo <- read.csv(paste0('../../misc/', cut, '/tsv/rowinfo.csv'),
                    row.names = 1)
dat_uncut <- read.table(paste0("../../data/", uncut, "/tsv/", uncut, ".tsv.gz"), 
                        sep="\t", 
                        header=T)
rowinfo_uncut <- read.csv(paste0('../../misc/', uncut, '/tsv/doublet_scores_DoubletDetection.csv'),
                    row.names = 2)

doublet <- read.csv(paste0('../../misc/', 
                           uncut, '/tsv/doublet_scores_DoubletDetection.csv'))
doublets <- c(doublet$Barcode[which(doublet$DoubletDetectionLabel==1)])
dat_uncut <- dat_uncut[!(row.names(dat_uncut)%in%doublets), ]
ampli_info <- read.table('../../misc/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
in_all <- intersect(colnames(dat_cut), intersect(colnames(dat_uncut), row.names(ampli_info)))
ampli_info <- ampli_info[][in_all, ]

fpr <- apply(dat_cut[row.names(rowinfo), row.names(ampli_info)[grepl("CpG.always.unmeth.B", ampli_info$Type.of.amplicon)]], 2, function(x){
  sum(x>0)/length(x)
})
fnr <- apply(dat_uncut[row.names(rowinfo_uncut[which(rowinfo_uncut$DoubletDetectionLabel==0), ]), row.names(ampli_info)[grepl("CpG.B.cell.diff", ampli_info$Type.of.amplicon)]], 2, function(x){
  sum(x==0)/length(x)
})
to_plot <- data.frame(Type=c(rep('FPR', length(fpr)),
                             rep('FNR', length(fnr))),
                     Value=c(fpr, fnr))
plot <- ggplot(to_plot, aes(x=Type, y=Value))+geom_boxplot(color='black', fill='gray80', size=.25, outlier.size=.25)+plot_theme+
  xlab('')+ylab('')
ggsave(file.path(plot.path, 'FPR_FNR_plot.pdf'), 
       plot,
       width=40,
       height=40,
       units='mm')
