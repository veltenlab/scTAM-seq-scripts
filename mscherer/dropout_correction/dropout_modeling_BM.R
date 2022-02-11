library(ggplot2)
library(rstan)
plot_path <- '/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure1/Sample11/correction/clusters/'
cut <- 'Sample11_70_percent_good_performance'
uncut <- 'Sample12_70_percent_good_performance'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=15),
                    axis.text=element_text(color='black',size=15),
                    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    legend.key=element_rect(color='black', fill=NA))
filtered.counts <- read.table(paste0("/users/lvelten/project/Methylome/analysis/missionbio/BM/", cut, "/tsv/", cut, ".barcode.cell.distribution.tsv"), 
                              sep="\t", 
                              header=T)
cell_metadata <- read.csv(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/', cut, '/tsv/rowinfo.csv'),
                          row.names = 1)
cell_metadata <- cell_metadata[]
all.amplicons <- read.table('/users/lvelten/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
# Only using the best-performing amplicons
#sel.amplicons <- row.names(amplicon.info)
# Using all amplicons
sel.amplicons <- row.names(all.amplicons[grepl('CpG.B.cell.diff', all.amplicons$Type.of.amplicon),])
non_cut <- read.table(paste0("/users/lvelten/project/Methylome/analysis/missionbio/BM/", uncut, "/tsv/", uncut, ".barcode.cell.distribution.tsv"), 
                      sep="\t", 
                      header=T)
non_cut_doublets <- read.csv(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/', 
                                    uncut, '/tsv/doublet_scores_DoubletDetection.csv'))
non_cut <- non_cut[non_cut_doublets[which(non_cut_doublets$DoubletDetectionLabel==0), 'Barcode'], ]
allelic_dropout <- apply(non_cut[, sel.amplicons], 2, function(x){
  sum(x==0)/length(x)
})
sel.amplicons <- sel.amplicons[allelic_dropout<=0.025&allelic_dropout>=0.01]
allelic_dropout <- allelic_dropout[sel.amplicons]

selected <- ifelse(filtered.counts[row.names(cell_metadata), sel.amplicons]>0, 1, 0)

vals_naive <- apply(selected[cell_metadata$CellType_detailed%in%'naive B-cells', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})
vals_memory <- apply(selected[cell_metadata$CellType_detailed%in%'memory B-cells', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})
vals_S1 <- apply(selected[cell_metadata$CellType_detailed%in%'S1 cells', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})
vals_S2 <- apply(selected[cell_metadata$CellType_detailed%in%'S2 cells', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})
vals_S3S4 <- apply(selected[cell_metadata$CellType_detailed%in%'S3-S4 cells', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})

corrected_values <- sapply(sel.amplicons, function(ampli){
  data_naive <- list(n0=vals_naive[1, ampli], n1=vals_naive[2, ampli], p=allelic_dropout[ampli])
  sm_naive <- stan("basic.stan", data = data_naive)
  data_memory <- list(n0=vals_memory[1, ampli], n1=vals_memory[2, ampli], p=allelic_dropout[ampli])
  sm_memory <- stan("basic.stan", data = data_memory)
  data_S1 <- list(n0=vals_S1[1, ampli], n1=vals_S1[2, ampli], p=allelic_dropout[ampli])
  sm_S1 <- stan("basic.stan", data = data_S1)
  data_S2 <- list(n0=vals_S2[1, ampli], n1=vals_S2[2, ampli], p=allelic_dropout[ampli])
  sm_S2 <- stan("basic.stan", data = data_S2)
  data_S3S4 <- list(n0=vals_S3S4[1, ampli], n1=vals_S3S4[2, ampli], p=allelic_dropout[ampli])
  sm_S3S4 <- stan("basic.stan", data = data_S3S4)
  c(NaiveB=get_posterior_mean(sm_naive)['m', 'mean-all chains'],
    MemoryB=get_posterior_mean(sm_memory)['m', 'mean-all chains'],
    S1cells=get_posterior_mean(sm_S1)['m', 'mean-all chains'],
    S2cells=get_posterior_mean(sm_S2)['m', 'mean-all chains'],
    S3S4cells=get_posterior_mean(sm_S3S4)['m', 'mean-all chains'])
})
