############## F1_basic_statistics.R ##############
#' This file generates basic information about the experiment, including false positive rate, false
#' negative rate and dropout rates

library(ggplot2)
cut <- 'Sample5_80_percent'
uncut <- 'Sample2_default'
plot.path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/revision/'
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
dat_cut <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/", cut, "/tsv/", cut, ".barcode.cell.distribution.tsv"), 
                      sep="\t", 
                      header=T)
dat_uncut <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/", uncut, "/tsv/", uncut, ".barcode.cell.distribution.tsv"), 
                        sep="\t", 
                        header=T)
fnr <- apply(dat_uncut, 2, function(x){
  sum(x==0)/length(x)
})
jurkat_high <- as.character(read.csv('/users/mscherer/cluster/project/Methylome/analysis/amplicons/cell_lines/confident_CpGs_Jurkat_high.csv', header = F)$V1)
k562_high <- as.character(read.csv('/users/mscherer/cluster/project/Methylome/analysis/amplicons/cell_lines/confident_CpGs_K562_high.csv', header = F)$V1)
dat_cut <- ifelse(dat_cut>0, 1, 0)
clust <- read.table('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/Sample5_80_percent/tsv/cluster_assignment.tsv')
fpr.jurkat <- colMeans(dat_cut[row.names(subset(clust, subset=CellType=='Jurkat')), k562_high])
fpr.k562 <- colMeans(dat_cut[row.names(subset(clust, subset=CellType=='K562')), jurkat_high])
fpr <- c(fpr.jurkat,
         fpr.k562)
to_plot <- data.frame(Type=c(rep('FPR', length(fpr)),
                             rep('FNR', length(fnr))),
                     Value=c(fpr, fnr))
plot <- ggplot(to_plot, aes(x=Type, y=Value))+geom_boxplot(color='black', fill='gray80', size=.25, outlier.size=.25)+plot_theme+
  xlab('')+ylab('')
ggsave('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/revision/FPR_FNR_plot.pdf', 
       plot,
       width=40,
       height=40,
       units='mm')
