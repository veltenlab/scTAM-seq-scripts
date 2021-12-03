library(RnBeads)
genes <- rnb.get.annotation('genes', assembly='hg19')
cpgs <- rnb.get.annotation('CpG', assembly = 'hg19')
cpg.names <- row.names(rnb.annotation2data.frame(cpgs))
cpgs <- unlist(cpgs)
names(cpgs) <- cpg.names

sel.cpgs <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/differential/differential_CpGs_ Cluster2avsCluster2c .csv')
sel.cpgs <- cpgs[sel.cpgs$CpGID]

closest.genes <- c()
for(i in 1:length(sel.cpgs)){
  closest.genes <- c(closest.genes, names(genes)[which.min(distance(sel.cpgs[i], cpgs))])
}