#group cells into memory, naive or potential doublets (According to autoencoder)
.libPaths(c(.libPaths(), '/users/mscherer/conda/envs/rnbeads/lib/R/library/'))
require(umap)
require(Seurat)
require(ggplot2)
require(plyr)
require(reshape2)
library(pheatmap)
filtered.counts <- read.table("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/tsv/BCells_Sample7_70_percent_good_performance.barcode.cell.distribution_with_MCL.tsv", row.names = 1, header=T)
amplicon.info <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv", header=T, row.names = 1)

bulk.methylation <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt", header=T)
doublet_file <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/doublet_detection/Hha_dropout/doublets.csv', row.names = 1)

autoencoder.bottleneck <- read.csv("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/methylation_autoencoder/bottleneck_good.csv", row.names = 1)
autoencoder.output <- read.csv("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/methylation_autoencoder/mixture_prob_dca_with_MCL.csv", row.names = 1)
missing.amplicons <- read.csv("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/methylation_autoencoder/missing_amplicons_with_MCL.csv", row.names = 1)
filtered.counts <- filtered.counts[,!colnames(filtered.counts) %in% missing.amplicons$X0]
dimnames(autoencoder.output) <- dimnames(filtered.counts)

colnames(autoencoder.bottleneck) <- paste0("BOTTLE_", 1:ncol(autoencoder.bottleneck))
rownames(autoencoder.bottleneck) <- rownames(filtered.counts)

amplicon.info <- amplicon.info[colnames(filtered.counts),]

#Create a seurat object... lol 
methylome <- CreateSeuratObject(counts = t(filtered.counts), min.cells = 20, assay = "TAM" ) #only amplicons observed in at least 20 cells
methylome@assays$TAM@data <- as.matrix(t(autoencoder.output[,rownames(methylome)]))
methylome@reductions$bottleneck <- CreateDimReducObject(embeddings = as.matrix(autoencoder.bottleneck), assay = "TAM",key = "BOTTLE_")
methylome <- RunUMAP(methylome, dims = 1:ncol(autoencoder.bottleneck), reduction = "bottleneck")
methylome <- FindNeighbors(methylome, reduction  = "bottleneck", dims = 1:ncol(autoencoder.bottleneck))
methylome <- FindClusters(methylome, resolution = 0.2)

DimPlot(methylome)
methylome$log_count <- log10(methylome$nCount_TAM)
VlnPlot(methylome, c("log_count", "nFeature_TAM")) #no systematic bias in read depth between clusters, except that cl.2 has systematically less amplicons and also reads because there are less methylated sites
FeaturePlot(methylome, c("AMPL130853","AMPL130724","AMPL130910","AMPL131020"), slot = "counts")
labels <- c("0"="naive_1", "1" = "naive_2", "3" = "inbetween", "2" = "memory")
Idents(methylome) <- labels[as.character(Idents(methylome))]

#Similar data structure for uncut experiment
filtered.counts.uncut <- read.table("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample6_70_percent_good_performance/tsv/BCells_Sample6_70_percent_good_performance.barcode.cell.distribution.tsv", row.names = 1, header=T)
nreads <- apply(filtered.counts.uncut, 1,sum)

#How well does the model explain the data in each cell type? (avg likelihood per cell)
#for some cells this will not work
get_likelihood_bygroup <- function(amplicon, plot=T) {
  out <- data.frame(x = c(scale(log10(nreads)), scale(methylome$log_count)),
                    y = c(filtered.counts.uncut[,amplicon], GetAssayData(methylome,slot = "counts")[amplicon,]) ,
                    class = c(rep("uncut", nrow(filtered.counts.uncut)), rep("cut", ncol(methylome))),
                    celltype = c(rep("uncut", nrow(filtered.counts.uncut)), as.character(Idents(methylome)))
  )
  out$binary <- out$y > 0
  
  model <- glm(binary ~ x + celltype, data= out, family = "binomial")
  
  out$fit <- predict(model, out, type = "response")
  summarised <- ddply(out, "celltype", summarise, logLig.percell = sum( ifelse(binary, -log10(fit), -log10(1-fit)))/length(fit), mean.y = mean(y), fr.pos = mean(binary) )
  summarised$amplicon <- amplicon
  return(summarised)
}
allLik <- lapply(rownames(methylome), get_likelihood_bygroup)
allLik <- do.call(rbind, allLik)

allLikShort <- dcast(allLik,amplicon ~ celltype, value.var = "fr.pos")
allLikShort <- melt(allLikShort, id.vars = c("amplicon","uncut"))
amplicon.info$amplicon <- rownames(amplicon.info)

allLikShort <- merge(allLikShort, amplicon.info[,c("amplicon","GC_Content", "Type.of.amplicon")])
#allLikShort <- merge(allLikShort, bulk.methylation, all.x=T)
#compute dropout probability as a function of total read depth
FeatureScatter(methylome, "log_count", "AMPL130853", slot = "counts") + scale_y_log10() + geom_smooth()
FeatureScatter(methylome, "log_count", "AMPL130724", slot = "counts") + scale_y_log10()


#conclusion: False positive rate is near zero and fix, false negative rate can accurately be estimated
#fpr per amplicon
fpr <- subset(allLikShort, variable == "memory" , select = c("amplicon","uncut", "Type.of.amplicon"))
rownames(fpr) <- fpr$amplicon
amplicon.info <- amplicon.info[rownames(methylome),]
fpr <- fpr[rownames(amplicon.info),]
#select relevant amplicons
filtered.counts <- filtered.counts[,rownames(methylome)]
amplicon.info$Type.of.amplicon[is.na(amplicon.info$Type.of.amplicon)] <- "NA"
selected <- filtered.counts[,fpr$Type.of.amplicon == "CpG.B.cell.diff" & fpr$uncut > 0.9] > 0
selected <- apply(selected,2,as.numeric)
rownames(selected) <- rownames(filtered.counts)
non_hha_reads <- apply(filtered.counts[, row.names(amplicon.info[amplicon.info$Type.of.amplicon%in%'NonHhaI',])], 1, sum)
rowinfo <- data.frame(row.names = rownames(filtered.counts), 
                      ct = Idents(methylome), 
                      libsize = log10(methylome$nCount_TAM), 
                      non_hha_reads=non_hha_reads, 
                      nfeature = methylome$nFeature_TAM,
                      doublet = doublet_file[rownames(filtered.counts), 'x'])
pheatmap(selected, annotation_col = subset(fpr, select = "uncut"),annotation_row = rowinfo, clustering_distance_cols = "binary", clustering_distance_rows = "binary", show_rownames = F, show_colnames = F, cutree_rows = 4, clustering_method = "ward.D2")

# plot bulk vs. single-cell estimated level
means_clusters <- aggregate(selected, by=list(rowinfo$ct), mean)
bulk_table <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
means_bulk <- bulk_table[, c('NBC.mean', 'MBC.mean')]
colnames(means_bulk) <- c('naiveB', 'memoryB')
row.names(means_bulk) <- bulk_table$amplicon
means_bulk <- as.data.frame(t(means_bulk))
means_bulk$'Group.1' <- row.names(means_bulk)
to_plot <- rbind(means_clusters, means_bulk[, colnames(means_clusters)[colnames(means_clusters)%in%colnames(means_bulk)]])
row.names(to_plot) <- to_plot$Group.1
to_plot <- as.data.frame(t(to_plot[, -1]))
to_plot$Amplicon <- row.names(to_plot)
to_plot <- reshape2::melt(to_plot, id=c('Amplicon', 'naiveB', 'memoryB'))
colnames(to_plot)[4:5] <- c('Cluster', 'Methylation')
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
plot <- ggplot(to_plot, aes_string(x='naiveB', y='Methylation', color='Cluster'))+geom_point()+geom_smooth(method='lm',se=FALSE)+xlim(0, 1)+ylim(0,1)+
  facet_grid(Cluster~.)+theme
ggsave('/users/mscherer/cluster/project/Methylome/analysis/BCells/bulk_vs_SC/naive_B_Cells.png', plot)
plot <- ggplot(to_plot, aes_string(x='memoryB', y='Methylation', color='Cluster'))+geom_point()+geom_smooth(method='lm',se=FALSE)+xlim(0, 1)+ylim(0,1)+
  facet_grid(Cluster~.)+theme
ggsave('/users/mscherer/cluster/project/Methylome/analysis/BCells/bulk_vs_SC/memory_B_Cells.png', plot)

p.vals <- sapply(unique(rowinfo$ct), function(ct1){
  sapply(unique(rowinfo$ct), function(ct2){
    if(as.character(ct1)>=as.character(ct2)) return(NA)
    unlist(apply(filtered.counts, 2, function(x, c1, c2){
      sel_ct1 <- row.names(rowinfo)[rowinfo$ct%in%c1]
      sel_ct2 <- row.names(rowinfo)[rowinfo$ct%in%c2]
      sort(wilcox.test(x[sel_ct1], x[sel_ct2])$p.value)
    }, c1=ct1, c2=ct2))
  })
})
row.names(p.vals) <- colnames(p.vals) <- unique(rowinfo$ct)
sort(unlist(p.vals['naive_2', 'naive_1']))

clusters_heatmap <- cutree(hclust(dist(selected, 'binary'), method = 'ward.D2'), 4)
rowinfo$Cluster <- clusters_heatmap
pheatmap(selected, annotation_col = subset(fpr, select = "uncut"),annotation_row = rowinfo, clustering_distance_cols = "binary", clustering_distance_rows = "binary", show_rownames = F, show_colnames = F, cutree_rows = 4, clustering_method = "ward.D2")
p.vals <- sapply(unique(rowinfo$Cluster), function(ct1){
  sapply(unique(rowinfo$Cluster), function(ct2){
    if(as.character(ct1)>=as.character(ct2)) return(NA)
    unlist(apply(filtered.counts, 2, function(x, c1, c2){
      sel_ct1 <- row.names(rowinfo)[rowinfo$Cluster%in%c1]
      sel_ct2 <- row.names(rowinfo)[rowinfo$Cluster%in%c2]
      sort(wilcox.test(x[sel_ct1], x[sel_ct2])$p.value)
    }, c1=ct1, c2=ct2))
  })
})
row.names(p.vals) <- colnames(p.vals) <- unique(rowinfo$Cluster)
meth_diff <- sapply(unique(rowinfo$Cluster), function(ct1){
  sapply(unique(rowinfo$Cluster), function(ct2){
    if(as.character(ct1)>=as.character(ct2)) return(NA)
    unlist(apply(filtered.counts, 2, function(x, c1, c2){
      sel_ct1 <- row.names(rowinfo)[rowinfo$Cluster%in%c1]
      sel_ct2 <- row.names(rowinfo)[rowinfo$Cluster%in%c2]
      meth_first <- sum(x[sel_ct1]>0)/length(sel_ct1)
      meth_second <- sum(x[sel_ct2]>0)/length(sel_ct2)
      return(c(meth_first, meth_second))
    }, c1=ct1, c2=ct2))
  })
})
row.names(meth_diff) <- colnames(meth_diff) <- unique(rowinfo$Cluster)
load('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/meth.data.numeric.Rdata')

bcell_clusters <- unlist(p.vals[4,1])
bcell_clusters <- sort(bcell_clusters)[1:20]
bcell_clusters_diff <- t(as.data.frame(meth_diff[4,1]))
colnames(bcell_clusters_diff) <- c('MeanCluster1', 'MeanCluster4')
cell_assignment <- read.table('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/MBC_assignment.txt')
meth.data.numeric <- meth.data.numeric[,as.character(cell_assignment$V5)]
more_info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
row.names(more_info) <- more_info$amplicon
joint_names <- intersect(row.names(more_info), names(bcell_clusters))
bcell_clusters <- bcell_clusters[joint_names]
bcell_clusters_diff <- bcell_clusters_diff[joint_names,]
names(bcell_clusters) <- more_info[names(bcell_clusters), 'background.cpgs']
row.names(bcell_clusters_diff) <- more_info[row.names(bcell_clusters_diff), 'background.cpgs']
bcell_clusters_diff <- bcell_clusters_diff[!is.na(row.names(bcell_clusters_diff)),]
mean_classes <- aggregate(t(meth.data.numeric[row.names(bcell_clusters_diff), ]), by=list(cell_assignment$V2), mean)
bcell_clusters_diff <- t(bcell_clusters_diff)
bcell_clusters_diff <- data.frame(Group.1=row.names(bcell_clusters_diff), bcell_clusters_diff)
colnames(bcell_clusters_diff)[-1] <- paste0(colnames(bcell_clusters_diff)[-1], '/', 
                                            more_info[match(colnames(bcell_clusters_diff)[-1], more_info$background.cpgs), 'amplicon'])
colnames(mean_classes)[-1] <- paste0(colnames(mean_classes)[-1], '/', 
                                            more_info[match(colnames(mean_classes)[-1], more_info$background.cpgs), 'amplicon'])
to_plot <- data.frame(Type=c('Bulk', 'Bulk', 'SingleCell', 'SingleCell'),
                      rbind(mean_classes, 
                        bcell_clusters_diff))
to_plot <- reshape2::melt(to_plot, id=c('Group.1', 'Type'))
colnames(to_plot)[3:4] <- c('CpGID', 'Methylation')
#to_plot$variable <- factor(to_plot$variable, levels=names(bcell_clusters)[order(bcell_clusters)])

library(viridis)
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
plot <- ggplot(to_plot,aes(x=Group.1, y=CpGID, fill=Methylation))+geom_tile()+geom_text(aes(label=format(Methylation, digits = 2)), color='white')+
  scale_fill_viridis()+facet_wrap(Type~., scales='free_x', strip.position = "top")+theme+xlab('')

memory_selected <- filtered.counts[row.names(rowinfo[rowinfo$Cluster%in%c(1),]),]
#memory_selected <- memory_selected[,  joint_names]
memory_selected <- as.data.frame(apply(memory_selected, 1, function(x){
  ifelse(x>0, 1, 0)
}))
pheatmap(t(memory_selected), 
         annotation_col = subset(fpr, select = "uncut"),
         annotation_row = rowinfo, 
         clustering_distance_cols = "binary", 
         clustering_distance_rows = "binary", 
         show_rownames = F, 
         show_colnames = TRUE, 
         clustering_method = "ward.D2")

clusters <- cutree(hclust(dist(t(memory_selected), 'binary'), method='ward.D2'),2)
sel_ct1 <- names(which(clusters==1))
sel_ct2 <- names(which(clusters==2))
res <- apply(filtered.counts, 2, function(x){
  meth_first <- sum(x[sel_ct1]>5)/length(sel_ct1)
  meth_second <- sum(x[sel_ct2]>5)/length(sel_ct2)
  return(c(wilcox.test(x[sel_ct1], x[sel_ct2])$p.value, meth_first, meth_second))
})
res <- t(as.data.frame(res))
colnames(res) <- c('P-value', 'MeanCluster1', 'MeanCluster2')
res <- res[order(res[,'P-value']),]
total_reads_c1 <- apply(filtered.counts[sel_ct1, row.names(amplicon.info[amplicon.info$Type.of.amplicon%in%'NonHhaI',])], 2, sum)
total_reads_c2 <- apply(filtered.counts[sel_ct2, row.names(amplicon.info[amplicon.info$Type.of.amplicon%in%'NonHhaI',])], 2, sum)
to_plot <- data.frame(Cluster1=total_reads_c1, Cluster2=total_reads_c2)
to_plot <- reshape2::melt(to_plot)
plot <- ggplot(to_plot, aes(x=variable,y=value))+geom_boxplot()
