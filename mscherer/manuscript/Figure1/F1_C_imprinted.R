############## F1_C_imprinted.R ##############
#' With this file, we plot the percentage of droput (square root, because two alleles exist) of the undigested
#' in comparison to the digested sample for the imprinted amplicons. We also compute whether the square fit
#' is better than a linear fit.
#' 
library(ggplot2)
cut <- 'BCells_Sample7_70_percent_good_performance'
uncut <- 'BCells_Sample6_70_percent_good_performance'
plot.path <- '/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure1/'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=15),
                    axis.text=element_text(color='black',size=15),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.key=element_rect(color='black', fill=NA))
colors_amplicons <- c("Mutation only"="#e5c494",
                      "Differential CpG Bcells"="#fc8d62",
                      "MCL-specific CpG"="#ffd92f",
                      "Differential CpG Bcells AND MCL-specific CpG"="#fc8d62",
                      "Always unmethylated CpG in Bcells"="#a6d854",
                      "Always methylated CpG in Bcells"="#66c2a5",
                      "Imprinted CpG"="#e78ac3",
                      "Imprinted CpG multiple"="#e78ac3",
                      "No HhaI cutsite"="#b3b3b3")
dat_cut <- read.table(paste0("/users/lvelten/project/Methylome/analysis/missionbio/tapestri/", cut, "/tsv/", cut, ".barcode.cell.distribution_with_MCL.tsv"), 
                      sep="\t", 
                      header=T)
dat_uncut <- read.table(paste0("/users/lvelten/project/Methylome/analysis/missionbio/tapestri/", uncut, "/tsv/", uncut, ".barcode.cell.distribution.tsv"), 
                        sep="\t", 
                        header=T)
ampli_info <- read.table('/users/lvelten/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
sel_amplis <- row.names(ampli_info)[grepl('CpG.imprinted', ampli_info$Type.of.amplicon)]
dropout_cut <- apply(dat_cut[, sel_amplis], 2, function(x){
  sum(x==0)/length(x)
})
dropout_uncut <- apply(dat_uncut[, sel_amplis], 2, function(x){
  sum(x==0)/length(x)
})
to_plot <- data.frame(Cut=dropout_cut,
                      Uncut=dropout_uncut)
to_plot$Type <- 'Observed'
square_fit <- lm(Cut~0+Uncut+I(sqrt(Uncut)), data=to_plot)
linear_fit <- lm(Cut~0+Uncut, data=to_plot)
anov <- anova(linear_fit, square_fit)
anov_p <- anov$'Pr(>F)'[2]

plot <- ggplot(to_plot, aes(x=sqrt(Uncut), y=Cut))+
  geom_point(color=unname(colors_amplicons['Imprinted CpG']))+geom_abline(slope=1, intercept=0)+
  labs(x='Undigested: sqrt(dropout rate)', y='Digested: dropout rate', col='Amplicon Type')+
  xlim(0,1)+ylim(0,1)+
  annotate(x=0.3, y=0.75, geom = 'text', label=paste("ANOVA p-value sqrt- vs. linear fit: ", format(anov_p, digits = 3)))+
  plot_theme
ggsave(file.path(plot.path, 'F1_C_imprinted_original.pdf'), plot)

square_fit <- lm(Cut~0+Uncut+I(Uncut^2), data=to_plot)
linear_fit <- lm(Cut~0+Uncut, data=to_plot)
anov <- anova(linear_fit, square_fit)
anov_p <- anov$'Pr(>F)'[2]

square <- data.frame(Uncut=seq(0,1,by=0.01), Cut=seq(0,1,by=0.01))
square$Type <- 'Theoretical'
to_plot <- as.data.frame(rbind(to_plot, square))
to_plot <- to_plot[order(to_plot$Type, decreasing = TRUE),]
row.names(to_plot) <- NULL
plot <- ggplot(to_plot, aes(x=1-(Uncut^2), y=1-Cut, color=Type))+
  geom_point()+geom_abline(slope=1, intercept=0)+
  labs(x='Undigested: 1-(dropout rate)²', y='Digested: 1-(dropout rate)', col='Amplicon Type')+
  xlim(0,1)+ylim(0,1)+
  annotate(x=0.3, y=0.75, geom = 'text', label=paste("ANOVA p-value squared- vs. linear fit: ", format(anov_p, digits = 3)))+
  scale_color_manual(values=c('Observed'=unname(colors_amplicons['Imprinted CpG']), 'Theoretical'='#b3b3b3'))+
  plot_theme
ggsave(file.path(plot.path, 'F1_C_imprinted.pdf'), plot)
