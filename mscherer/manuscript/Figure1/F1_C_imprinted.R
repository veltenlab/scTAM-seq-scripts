############## F1_C_imprinted.R ##############
#' With this file, we plot the percentage of droput (square root, because two alleles exist) of the undigested
#' in comparison to the digested sample for the imprinted amplicons. We also compute whether the square fit
#' is better than a linear fit.
#' 
library(ggplot2)
cut <- 'Sample8_70_percent_good_performance'
uncut <- 'Sample12_70_percent_good_performance'
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
rowinfo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', cut, '/tsv/rowinfo.csv'),
                    row.names = 1)
dat_cut <- dat_cut[row.names(rowinfo), ]
dat_uncut <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/BM/", uncut, "/tsv/", uncut, ".barcode.cell.distribution.tsv"), 
                        sep="\t", 
                        header=T)
doublet <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/BM/', 
                           uncut, '/tsv/doublet_scores_DoubletDetection.csv'))
doublets <- c(doublet$Barcode[which(doublet$DoubletDetectionLabel==1)])
dat_uncut <- dat_uncut[!(row.names(dat_uncut)%in%doublets), ]
ampli_info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
sel_amplis <- row.names(ampli_info)[grepl('CpG.imprinted', ampli_info$Type.of.amplicon)]
dropout_cut <- apply(dat_cut[, sel_amplis], 2, function(x){
  sum(x>0)/length(x)
})
dropout_uncut <- apply(dat_uncut[, sel_amplis], 2, function(x){
  sum(x>0)/length(x)
})
to_plot <- data.frame(Cut=dropout_cut,
                      Uncut=dropout_uncut)
to_plot$Type <- 'Observed'
square_fit <- lm(Cut~0+Uncut+I(sqrt(Uncut)), data=to_plot)
linear_fit <- lm(Cut~0+Uncut, data=to_plot)
anov <- anova(linear_fit, square_fit)
anov_p <- anov$'Pr(>F)'[2]

plot <- ggplot(to_plot, aes(x=Uncut, y=Cut))+
  geom_point(color=unname(colors_amplicons['Imprinted CpG']))+geom_abline(slope=1, intercept=0)+
  labs(x='Fraction of cells with reads undigested sample', y='Fraction of cells with reads digested sample', col='Amplicon Type')+
  xlim(0,1)+ylim(0,1)+
  annotate(x=0.3, y=0.75, geom = 'text', label=paste("ANOVA p-value sqrt- vs. linear fit: ", format(anov_p, digits = 3)))+
  plot_theme
ggsave(file.path(plot.path, 'F1_C_imprinted_original.png'), plot, width=65, height=45, units='mm')

square_fit <- lm(Cut~0+Uncut+I(sqrt(Uncut)), data=to_plot)
linear_fit <- lm(Cut~0+Uncut, data=to_plot)
anov <- anova(linear_fit, square_fit)
anov_p <- anov$'Pr(>F)'[2]

square <- data.frame(Uncut=sqrt(seq(0,1,by=0.001)), Cut=seq(0,1,by=0.001))
#square$Type <- 'Sqrt'
#to_plot <- as.data.frame(rbind(to_plot, square))
#to_plot <- to_plot[order(to_plot$Type, decreasing = TRUE),]
row.names(to_plot) <- NULL
plot <- ggplot()+
  geom_point(data = square, aes(x=Uncut, y=Cut), size=.2, color='black')+
  geom_point(data = to_plot, aes(x=Uncut, y=Cut), color=colors_amplicons["Imprinted CpG"])+
  geom_abline(slope=1, intercept=0)+
  labs(x='Fraction of cells with reads undigested sample', y='Fraction of cells with reads digested sample', col='Amplicon Type')+
  xlim(0,1)+ylim(0,1)+
  annotate(x=0.45, y=1.0, geom = 'text', label=paste("ANOVA p-value sqrt- vs. linear fit: ", format(anov_p, digits = 3)), size=2)+
  scale_color_manual(values=c('Observed'=unname(colors_amplicons['Imprinted CpG']), 'Sqrt'='black'))+
  plot_theme
ggsave(file.path(plot.path, 'F1_C_imprinted.pdf'), plot, width=65, height=45, units='mm')
