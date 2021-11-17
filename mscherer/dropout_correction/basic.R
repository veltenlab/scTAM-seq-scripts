library(rstan)
cut <- 'BCells_Sample7_70_percent_good_performance'
uncut <- 'BCells_Sample6_70_percent_good_performance'
dat_cut <- read.table(paste0("/users/lvelten/project/Methylome/analysis/missionbio/tapestri/", cut, "/tsv/", cut, ".barcode.cell.distribution_with_MCL.tsv"), 
                      sep="\t", 
                      header=T)
rowinfo <- read.csv('/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure1/rowinfo_reclustering_doublet.csv',
                    row.names = 1)
dat_cut <- dat_cut[row.names(rowinfo), ]
dat_uncut <- read.table(paste0("/users/lvelten/project/Methylome/analysis/missionbio/tapestri/", uncut, "/tsv/", uncut, ".barcode.cell.distribution.tsv"), 
                        sep="\t", 
                        header=T)
doublet <- read.csv('/users/lvelten/project/Methylome/analysis/missionbio/tapestri/BCells_Sample6_70_percent_good_performance/tsv/doublet_scores_DoubletDetection.csv')
doublets <- c(doublet$Barcode[doublet$DoubletDetectionLabel==1])
dat_uncut <- dat_uncut[!(row.names(dat_uncut)%in%doublets), ]
ampli_info <- read.table('/users/lvelten/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
sel_amplicons <- row.names(ampli_info)[ampli_info$Type.of.amplicon%in%'CpG.always.meth.B']
dropout_uncut <- apply(dat_uncut, 2, function(x){
    sum(x==0)/length(x)
})
vals <- apply(dat_cut, 2, function(x){
    c(sum(x==0), sum(x>0))
})
corrected_values <- sapply(sel_amplicons, function(sel_amplicon){
    data <- list(n0=vals[1, sel_amplicon],n1=vals[2, sel_amplicon], p=dropout_uncut[sel_amplicon])
    sm <- stan("basic.stan",data = data)
    get_posterior_mean(sm)['m', 'mean-all chains']
})
uncorrected_values <- apply(dat_cut[, sel_amplicons], 2, function(x){
    sum(x>0)/length(x)
})

