library(ggplot2)
sample <- 'Sample3_default'
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
scrub <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores.csv'))
row.names(scrub) <- scrub$Barcode
scrub <- scrub[row.names(clust.file),]
table(scrub$PredictedDoublet,clust.file$CellType)
to.plot <- data.frame(CellType=clust.file$CellType,DoubletScore=scrub$DoubletScore)
plot <- ggplot(to.plot,aes(x=CellType,y=DoubletScore))+geom_boxplot()+theme_bw()
ggsave(paste0('/users/mscherer/cluster/project/Methylome/analysis/doublet_detection/',sample,'label_vs_doublet_score.png'))
