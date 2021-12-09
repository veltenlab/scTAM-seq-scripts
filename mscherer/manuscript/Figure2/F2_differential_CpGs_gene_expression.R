################### F2_differential_CpGs_gene_expression.R ################################
#' With this script, we use the differential CpGs defined before, located the closest gene
#' and investigate the gene expression in cs- and ncs-memory Bcells in Sergio's data
#' (10.1038/s41590-021-01059-0) for those CpGs that are within the gene body (+-5kb).

library(RnBeads)
library(Seurat)
genes <- unlist(rnb.get.annotation('genes', assembly='hg19'))
genes <- resize(genes, width = width(genes)+(5000), fix = "center")
cpgs <- rnb.get.annotation('probes450', assembly = 'hg19')
cpg.names <- row.names(rnb.annotation2data.frame(cpgs))
cpgs <- unlist(cpgs)
names(cpgs) <- cpg.names

plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=10),
                    axis.text=element_text(color='black',size=8),
                    axis.text.x=element_text(angle=90, hjust=1),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')

sel.cpgs <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/differential/differential_CpGs_Cluster2avsCluster2b.csv')
sel.cpgs <- cpgs[sel.cpgs$CpGID]

closest.genes <- c()
closest.distance <- c()
for(i in 1:length(sel.cpgs)){
  disti <- distance(sel.cpgs[i], genes)
  ns <- names(closest.distance)
  closest.distance <- c(closest.distance, min(disti, na.rm=TRUE))
  names(closest.distance) <- c(ns, values(genes)[which.min(disti), 'symbol'])
}
sel_genes <- na.omit(names(closest.distance)[which(closest.distance==0)])
seurat.obj <- readRDS('/users/mscherer/cluster/project/AML/gene_expression/data/Sergio_figshare/WTA_projected.rds')
gene_exp_data <- GetAssayData(seurat.obj, slot='data')
idents <- Idents(seurat.obj)
sel_data <- gene_exp_data[sel_genes[sel_genes%in%row.names(gene_exp_data)], idents%in%c('Class-switched memory B cells', 'Unswitched memory B cells')]
to_plot <- as.data.frame(t(as.data.frame(sel_data)))
idents <- idents[idents%in%c('Class-switched memory B cells', 'Unswitched memory B cells')]
wilcox.p <- apply(to_plot, 2, function(x){
  c1 <- x[idents%in%'Class-switched memory B cells']
  c2 <- x[idents%in%'Unswitched memory B cells']
  c1 <- c1[c1>0]
  c2 <- c2[c2>0]
  if(length(c1)<2|length(c2)<2) return(NA)
  wilcox.test(c1, c2)$p.value
})
map <- c('Class-switched memory B cells'='class-switched memory B-cell',
         'Unswitched memory B cells'='non-switch memory B-cell')
to_plot$CellType <- map[as.character(idents)]
to_plot <- reshape2::melt(to_plot, id='CellType')
colnames(to_plot)[2:3] <- c('Gene', 'NormalizedExpression')
#to_plot <- to_plot[to_plot$NormalizedExpression>0, ]
plot <- ggplot(to_plot, aes(x=CellType, y=NormalizedExpression))+geom_violin()+facet_wrap(Gene~.)+plot_theme

FeaturePlot(seurat.obj, features=as.character(unique(to_plot$Gene)), reduction = 'Projected')
