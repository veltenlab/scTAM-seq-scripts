library(ggplot2)
library(rstan)
plot_path <- '/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure1/'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=15),
                    axis.text=element_text(color='black',size=15),
                    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    legend.key=element_rect(color='black', fill=NA))
filtered.counts <- read.table("/users/lvelten/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/tsv/BCells_Sample7_70_percent_good_performance.barcode.cell.distribution_with_MCL.tsv", row.names = 1, header=T)
cell_metadata <- read.csv('/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure1/rowinfo_reclustering_doublet.csv',
                          row.names=1)
amplicon.info <- read.csv("/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure1/selected_amplicons.csv", row.names = 1)
all.amplicons <- read.table('/users/lvelten/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
# Only using the best-performing amplicons
#sel.amplicons <- row.names(amplicon.info)
# Using all amplicons
sel.amplicons <- row.names(all.amplicons[grepl('CpG.B.cell.diff', all.amplicons$Type.of.amplicon),])
non_cut <- read.table("/users/lvelten/project/Methylome/analysis/missionbio/tapestri/BCells_Sample6_70_percent_good_performance/tsv/BCells_Sample6_70_percent_good_performance.barcode.cell.distribution.tsv", row.names = 1, header=T)
non_cut_doublets <- na.omit(read.csv("/users/lvelten/project/Methylome/analysis/missionbio/tapestri/BCells_Sample6_70_percent_good_performance/tsv/doublet_scores_DoubletDetection.csv"))
non_cut <- non_cut[non_cut_doublets[non_cut_doublets$DoubletDetectionLabel==0, 'Barcode'], ]
allelic_dropout <- apply(non_cut[, sel.amplicons], 2, function(x){
  sum(x==0)/length(x)
})
sel.amplicons <- sel.amplicons[allelic_dropout>=0.75]
allelic_dropout <- allelic_dropout[sel.amplicons]

selected <- ifelse(filtered.counts[row.names(cell_metadata), sel.amplicons]>0, 1, 0)

vals_c1 <- apply(selected[cell_metadata$CellType_reclustering%in%'Cluster1', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})
vals_c2a <- apply(selected[cell_metadata$CellType_reclustering%in%'Cluster2a', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})
vals_c2b <- apply(selected[cell_metadata$CellType_reclustering%in%'Cluster2b', ], 2, function(ampli){
  c(sum(ampli==0),
    sum(ampli>0),
    sum(ampli>0)/length(ampli))
})
vals_c2c <- apply(selected[cell_metadata$CellType_reclustering%in%'Cluster2c', ], 2, function(ampli){
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
  data_c2c <- list(n0=vals_c2c[1, ampli], n1=vals_c2c[2, ampli], p=allelic_dropout[ampli])
  sm_c2c <- stan("basic.stan", data = data_c2c)
  c(Cluster1=get_posterior_mean(sm_c1)['m', 'mean-all chains'],
    Cluster2a=get_posterior_mean(sm_c2a)['m', 'mean-all chains'],
    Cluster2b=get_posterior_mean(sm_c2b)['m', 'mean-all chains'],
    Cluster2c=get_posterior_mean(sm_c2c)['m', 'mean-all chains'])
})
