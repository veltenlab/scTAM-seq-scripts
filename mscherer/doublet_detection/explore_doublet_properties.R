library(ggplot2)
sample <- 'Sample5_80_percent'
dat.s3 <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'))
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
prot.counts <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample, '-protein-counts.tsv'),sep='\t')
ampli.info <- read.csv('/users/mscherer/cluster/project/Methylome/infos/cell_lines/analysis_Agostina_Apr21/Summary_methylation_values_amplicons_K562_Jurkat_v2.csv')
dat.s3 <- dat.s3[as.character(ampli.info$AmpID_design_1898)]
total.counts <- apply(dat.s3[,ampli.info$Type%in%"Aci"],1,sum)
to.plot <- data.frame(TotalCount=total.counts[row.names(clust.file)],Label=clust.file$CellType)
plot <- ggplot(to.plot,aes(x=Label,y=TotalCount))+geom_boxplot()+theme_bw()
ggsave(paste0('/users/mscherer/cluster/project/Methylome/analysis/doublet_detection/',sample,'label_vs_reads.png'))
total.protein.counts <- apply(prot.counts,1,sum)
to.plot <- data.frame(TotalCount=total.protein.counts[gsub('-1','',row.names(clust.file))],Label=clust.file$CellType)
plot <- ggplot(to.plot,aes(x=Label,y=log10(TotalCount+1)))+geom_boxplot()+theme_bw()
ggsave(paste0('/users/mscherer/cluster/project/Methylome/analysis/doublet_detection/',sample,'label_vs_protein_reads.png'))
t.test(to.plot$TotalCount[to.plot$Label=='Jurkat'],to.plot$TotalCount[to.plot$Label=='Mixed'])
t.test(to.plot$TotalCount[to.plot$Label=='K562'],to.plot$TotalCount[to.plot$Label=='Mixed'])
sum.mixed <- sum(to.plot$Label=='Mixed')
sum.jurkat <- sum(to.plot$Label=='Jurkat')
sum.k562 <- sum(to.plot$Label=='K562')
bartlett.test(c(to.plot$TotalCount[to.plot$Label=='Jurkat'],to.plot$TotalCount[to.plot$Label=='Mixed']),
              g=c(rep('Jurkat',sum.jurkat),rep('Mixed',sum.mixed)))
bartlett.test(c(to.plot$TotalCount[to.plot$Label=='K562'],to.plot$TotalCount[to.plot$Label=='Mixed']),
              g=c(rep('K562',sum.k562),rep('Mixed',sum.mixed)))
