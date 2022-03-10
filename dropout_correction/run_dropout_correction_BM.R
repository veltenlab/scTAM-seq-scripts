############################ run_dropout_correction.R ############################ 
#' This file is used to estimate pseudo-bulk DNA-methylation values for cell type
#' cluster (in the cut sample) accounting for the dropout rate in an uncut experiment

library(ggplot2)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
plot_path <- '~'
cut <- 'GSM5935921_BM_HhaI'
uncut <- 'GSM5935923_BM_undigested'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=15),
                    axis.text=element_text(color='black',size=15),
                    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    legend.key=element_rect(color='black', fill=NA))
filtered.counts <- read.table(paste0('../data/', cut, ".tsv.gz"), 
                              sep="\t", 
                              header=T)
cell_metadata <- read.csv(paste0('../misc/', cut, '/tsv/rowinfo.csv'),
                          row.names = 1)
cell_metadata <- cell_metadata[]
all.amplicons <- read.table('../misc/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
sel.amplicons <- row.names(all.amplicons[grepl('CpG.B.cell.diff', all.amplicons$Type.of.amplicon),])
non_cut <- read.table(paste0('../data/', uncut, ".tsv.gz"), 
                      sep="\t", 
                      header=T)
non_cut_doublets <- read.csv(paste0('../misc/', uncut, '/tsv/doublet_scores_DoubletDetection.csv'))
non_cut <- non_cut[non_cut_doublets[which(non_cut_doublets$DoubletDetectionLabel==0), 'Barcode'], ]
allelic_dropout <- apply(non_cut[, sel.amplicons], 2, function(x){
  sum(x==0)/length(x)
})

selected <- ifelse(filtered.counts[row.names(cell_metadata), sel.amplicons]>0, 1, 0)

vals_c1 <- apply(selected[cell_metadata$CellType_detailed%in%'S1 cells', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})
vals_c2 <- apply(selected[cell_metadata$CellType_detailed%in%'S2 cells', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})
vals_c3 <- apply(selected[cell_metadata$CellType_detailed%in%'S3-S4 cells', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})
vals_c4 <- apply(selected[cell_metadata$CellType_detailed%in%'naive B-cells', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})

vals_c5 <- apply(selected[cell_metadata$CellType_detailed%in%'memory B-cells', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})


corrected_values <- sapply(sel.amplicons, function(ampli){
  data_c1 <- list(n0=vals_c1[1, ampli], n1=vals_c1[2, ampli], p=allelic_dropout[ampli])
  sm_c1 <- stan("basic.stan", data = data_c1)
  data_c2 <- list(n0=vals_c2[1, ampli], n1=vals_c2[2, ampli], p=allelic_dropout[ampli])
  sm_c2 <- stan("basic.stan", data = data_c2)
  data_c3 <- list(n0=vals_c3[1, ampli], n1=vals_c3[2, ampli], p=allelic_dropout[ampli])
  sm_c3 <- stan("basic.stan", data = data_c3)
  data_c4 <- list(n0=vals_c4[1, ampli], n1=vals_c4[2, ampli], p=allelic_dropout[ampli])
  sm_c4 <- stan("basic.stan", data = data_c4)
  data_c5 <- list(n0=vals_c5[1, ampli], n1=vals_c5[2, ampli], p=allelic_dropout[ampli])
  sm_c5 <- stan("basic.stan", data = data_c5)
  c('S1 cells'=get_posterior_mean(sm_c1)['m', 'mean-all chains'],
    'S2 cells'=get_posterior_mean(sm_c2)['m', 'mean-all chains'],
    'S3-S4 cells'=get_posterior_mean(sm_c3)['m', 'mean-all chains'],
    'naive B-cells'=get_posterior_mean(sm_c4)['m', 'mean-all chains'],
    'memory B-cells'=get_posterior_mean(sm_c5)['m', 'mean-all chains'])
})
if(!dir.exists(file.path('../dropout_modeling'))){
  dir.create(file.path('../dropout_modeling'))
}
if(!dir.exists(file.path('../dropout_modeling', cut))){
  dir.create(file.path('../dropout_modeling', cut))
}
write.csv(t(corrected_values), file.path('../dropout_modeling', cut, 'corrected_values.csv'))