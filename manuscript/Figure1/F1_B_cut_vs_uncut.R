############## F1_B_cut_vs_uncut.R ##############
#' This file generates the plot comparing the dropout rate per amplicon (type) in the cut versus the uncut
#' sample. Dropout is defined as those cells having 0 reads at this particular amplicon

library(ggplot2)
cut <- 'GSM5935921_BM_HhaI'
uncut <- 'GSM5935923_BM_undigested'
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
dat_cut <- read.table(paste0('../../data/', cut, ".tsv.gz"), 
                  sep="\t", 
                  header=T)
rowinfo <- read.csv(paste0('../../misc/', cut, '/tsv/rowinfo.csv'), row.names=1)
non_doublets <- row.names(rowinfo)[which(rowinfo$DoubletDetectionLabel==0)]
dat_cut <- dat_cut[non_doublets, ]
dat_uncut <- read.table(paste0('../../data/', uncut, ".tsv.gz"), 
                      sep="\t", 
                      header=T)
doublet <- read.csv(paste0('../../misc/', uncut, '/tsv/doublet_scores_DoubletDetection.csv'))
non_doublets <- doublet$Barcode[which(doublet$DoubletDetectionLabel==0)]
dat_uncut <- dat_uncut[non_doublets, ]
ampli_info <- read.table('../../misc/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
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
cors <- c()
for(i in 1:nrow(to_plot)){
  ty <- to_plot[i, 'Type']
  cori <- cor(to_plot$Cut[to_plot$Type%in%ty], to_plot$Uncut[to_plot$Type%in%ty])
  cors <- c(cors, cori)
}
to_plot$Correlation <- cors

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
  geom_text(aes(label=paste('r:', format(Correlation, digits=2))), y=.95, x=0.2, geom='text', check_overlap=TRUE, color='black', size=2, fontface='bold')+
  labs(x='Fraction of cells with reads in undigested sample', y='Fraction of cells with reads in digested sample', col='Amplicon Type')+
  plot_theme+scale_color_manual(values=colors_amplicons)
ggsave(file.path(plot.path, 'F1_B_cut_vs_uncut_selected.pdf'), plot, width=120, height=38, units='mm')

selected_amplicons <- names(which(dropout_uncut>0.75 & ampli_info$Type.of.amplicon%in%'CpG.B.cell.diff'))
selected_amplicons <- ampli_info[selected_amplicons, ]                          
