library(scran)
sample <- 'Sample3_80_percent'
dat <- read.table(paste0('/users/lvelten/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'))
res <- doubletCells(t(dat))
write.csv(res,paste0('/users/lvelten/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores_DoubletCells.csv'))
