############################### Supplement_Bcell_heatmap_doublets.R ######################
#' This file generates the heatmap after filtering for amplicons that work well in the 
#' undigested Sample. 

.libPaths(c(.libPaths(), '/users/mscherer/conda/envs/rnbeads/lib/R/library/'))
library(Seurat)
library(ggplot2)
library(plyr)
library(reshape2)
library(pheatmap)
library(viridis)
library(grid)
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Supplement/'
color_map <- list(CellType=c('naive'='#fcbd7e',
                 'memory1'='#fc6571',
                 'memory2'='#fc3262'),
               Bulk=c('Naive B-cell high'='#fcbd7e',
               'Memory B-cell high'='#fc4762'))
filtered.counts <- read.table("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/tsv/BCells_Sample7_70_percent_good_performance.barcode.cell.distribution_with_MCL.tsv", row.names = 1, header=T)
amplicon.info <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv", header=T, row.names = 1)

bulk.methylation <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt", header=T)
doublet_file <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/tsv/doublet_scores_DoubletDetection.csv', row.names = 2)

autoencoder.bottleneck <- read.csv("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/methylation_autoencoder/bottleneck_good.csv", row.names = 1)
autoencoder.output <- read.csv("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/methylation_autoencoder/mixture_prob_dca_with_MCL.csv", row.names = 1)
missing.amplicons <- read.csv("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/methylation_autoencoder/missing_amplicons_with_MCL.csv", row.names = 1)
filtered.counts <- filtered.counts[,!colnames(filtered.counts) %in% missing.amplicons$X0]
dimnames(autoencoder.output) <- dimnames(filtered.counts)

colnames(autoencoder.bottleneck) <- paste0("BOTTLE_", 1:ncol(autoencoder.bottleneck))
rownames(autoencoder.bottleneck) <- rownames(filtered.counts)

amplicon.info <- amplicon.info[colnames(filtered.counts),]

methylome <- CreateSeuratObject(counts = t(filtered.counts), min.cells = 20, assay = "TAM" ) #only amplicons observed in at least 20 cells
methylome@assays$TAM@data <- as.matrix(t(autoencoder.output[,rownames(methylome)]))
methylome@reductions$bottleneck <- CreateDimReducObject(embeddings = as.matrix(autoencoder.bottleneck), assay = "TAM",key = "BOTTLE_")
methylome <- RunUMAP(methylome, dims = 1:ncol(autoencoder.bottleneck), reduction = "bottleneck")
methylome <- FindNeighbors(methylome, reduction  = "bottleneck", dims = 1:ncol(autoencoder.bottleneck))
methylome <- FindClusters(methylome, resolution = 0.2)

methylome$log_count <- log10(methylome$nCount_TAM)
FeaturePlot(methylome, c("AMPL130853","AMPL130724","AMPL130910","AMPL131020"), slot = "counts")
labels <- c("0"="naive_1", "1" = "naive_2", "3" = "inbetween", "2" = "memory")
Idents(methylome) <- labels[as.character(Idents(methylome))]

filtered.counts.uncut <- read.table("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample6_70_percent_good_performance/tsv/BCells_Sample6_70_percent_good_performance.barcode.cell.distribution.tsv", row.names = 1, header=T)
nreads <- apply(filtered.counts.uncut, 1,sum)

## How well does the model explain the data in each cell type? (avg likelihood per cell)
## for some cells this will not work
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
                      NonCutAmplicons=log10(non_hha_reads), 
                      Nfeatures = methylome$nFeature_TAM,
                      Doublet = ifelse(doublet_file[rownames(filtered.counts), 'DoubletDetectionLabel']==0, 'Singlet', 'Doublet'))

row.names(bulk.methylation) <- bulk.methylation$amplicon
bulk.methylation <- bulk.methylation[row.names(fpr),]
fpr$Bulk <- ifelse(bulk.methylation$NBC.mean>bulk.methylation$MBC.mean, 'Naive B-cell high', 'Memory B-cell high')

clusters <- cutree(hclust(dist(selected,
                               'binary'),
                      method='ward.D2'),
                   3)
cluster_map <- c('1'='memory2',
                 '2'='naive',
                 '3'='memory1')
rowinfo$CellType <- cluster_map[as.character(clusters)]
ph <- pheatmap(selected[row.names(rowinfo), ], 
               annotation_col = subset(fpr, select = c("Bulk")),
               annotation_row = subset(rowinfo, select = c("CellType", "Doublet")), 
               clustering_distance_cols = "binary", 
               clustering_distance_rows = "binary", 
               show_rownames = F, 
               show_colnames = F, 
               cutree_rows = 3, 
               clustering_method = "ward.D2",
               color=rev(inferno(50)),
               annotation_colors = color_map,
               fontsize=15)
ph$gtable$grobs[[1]]$gp <- gpar(lwd = 2)
ph$gtable$grobs[[2]]$gp <- gpar(lwd = 2)
png(file.path(plot_path, 'Supplement_BCell_heatmap.png'),
    width=1000,
    height=1600)
grid.newpage()
grid.draw(ph$gtable)
dev.off()