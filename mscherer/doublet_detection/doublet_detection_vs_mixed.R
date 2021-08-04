library(ggplot2)
sample <- 'Sample5_80_percent'
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
ddoublet <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/doublet_detection/DoubletDetection_classes.csv'))
ddoublet <- data.frame(ddoublet,read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/doublet_detection/DoubletDetection_scores.csv')))
row.names(ddoublet) <- row.names(clust.file)
ddoublet <- ddoublet[,-c(1,3)]
colnames(ddoublet) <- c('PredictedDoublet','DoubletScore')
table(ddoublet$PredictedDoublet,clust.file$CellType)
to.plot <- data.frame(CellType=clust.file$CellType,DoubletScore=ddoublet$DoubletScore)
plot <- ggplot(to.plot,aes(x=CellType,y=DoubletScore))+geom_boxplot()+theme_bw()
ggsave(paste0('/users/mscherer/cluster/project/Methylome/analysis/doublet_detection/',sample,'label_vs_doublet_score.png'))
