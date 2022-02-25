############################ run_dropout_correction.R ############################ 
#' This file is used to estimate pseudo-bulk DNA-methylation values for cell type
#' cluster (in the cut sample) accounting for the dropout rate in an uncut experiment

library(ggplot2)
library(rstan)
plot_path <- '~'
cut <- 'Sample8_70_percent_good_performance'
uncut <- 'Sample12_70_percent_good_performance'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=15),
                    axis.text=element_text(color='black',size=15),
                    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    legend.key=element_rect(color='black', fill=NA))
filtered.counts <- read.table(paste0(cut, "/tsv/", cut, ".barcode.cell.distribution.tsv"), 
                              sep="\t", 
                              header=T)
cell_metadata <- read.csv(paste0(cut, '/tsv/rowinfo.csv'),
                          row.names = 1)
cell_metadata <- cell_metadata[]
all.amplicons <- read.table('../misc/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
sel.amplicons <- row.names(all.amplicons[grepl('CpG.B.cell.diff', all.amplicons$Type.of.amplicon),])
non_cut <- read.table(paste0(uncut, "/tsv/", uncut, ".barcode.cell.distribution.tsv"), 
                      sep="\t", 
                      header=T)
non_cut_doublets <- read.csv(paste0(uncut, '/tsv/doublet_scores_DoubletDetection.csv'))
non_cut <- non_cut[non_cut_doublets[which(non_cut_doublets$DoubletDetectionLabel==0), 'Barcode'], ]
allelic_dropout <- apply(non_cut[, sel.amplicons], 2, function(x){
  sum(x==0)/length(x)
})

selected <- ifelse(filtered.counts[row.names(cell_metadata), sel.amplicons]>0, 1, 0)

vals_c1 <- apply(selected[cell_metadata$Cluster%in%'Cluster1', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})
vals_c2a <- apply(selected[cell_metadata$Cluster%in%'Cluster2a', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})
vals_c2b <- apply(selected[cell_metadata$Cluster%in%'Cluster2b', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})

corrected_values <- sapply(sel.amplicons, function(ampli){
  data_c1 <- list(n0=vals_c1[1, ampli], n1=vals_c1[2, ampli], p=allelic_dropout[ampli])
  sm_c1 <- stan("basic.stan", data = data_c1)
  data_c2a <- list(n0=vals_c2a[1, ampli], n1=vals_c2a[2, ampli], p=allelic_dropout[ampli])
  sm_c2a <- stan("basic.stan", data = data_c2a)
  data_c2b <- list(n0=vals_c2b[1, ampli], n1=vals_c2b[2, ampli], p=allelic_dropout[ampli])
  sm_c2b <- stan("basic.stan", data = data_c2b)
  c(Cluster1=get_posterior_mean(sm_c1)['m', 'mean-all chains'],
    Cluster2a=get_posterior_mean(sm_c2a)['m', 'mean-all chains'],
    Cluster2b=get_posterior_mean(sm_c2b)['m', 'mean-all chains'])
})
write.csv(t(corrected_values), file.path('dropout_modeling', cut, 'corrected_values.csv'))