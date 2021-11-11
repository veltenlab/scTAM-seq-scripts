############################### F1_D_heatmap.R ###########################################
#' This file generates the heatmap after filtering for amplicons that work well in the 
#' undigested Sample. Furthermore, doublets are filtered according to the DoubletDetection
#' output.
.libPaths(c(.libPaths(), '/users/mscherer/conda/envs/rnbeads/lib/R/library/'))
library(Seurat)
library(ggplot2)
library(plyr)
library(reshape2)
library(pheatmap)
library(viridis)
library(grid)
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Supplement/'
color_map <- list(CellType=c('Jurkat'='#ff9289',
                 'K562'='#82b7ff'))
filtered.counts <- read.table("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/Sample5_80_percent/tsv/Sample5_80_percent.barcode.cell.distribution.tsv", row.names = 1, header=T)
amplicon.info <- read.table("/users/mscherer/cluster/project/Methylome/infos/amplicon_info/200804_Design_Mission_Bio_sites_Jurkat_K562.tsv", header=T, row.names = 10)

doublet_file <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/Sample5_80_percent/tsv/doublet_scores_DoubletDetection.csv', row.names = 2)

autoencoder.bottleneck <- read.csv("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/Sample5_80_percent/methylation_autoencoder/bottleneck.csv", row.names = 1)
autoencoder.output <- read.csv("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/Sample5_80_percent/methylation_autoencoder/mixture_prob_dca.csv", row.names = 1)
dimnames(autoencoder.output) <- dimnames(filtered.counts)

colnames(autoencoder.bottleneck) <- paste0("BOTTLE_", 1:ncol(autoencoder.bottleneck))
rownames(autoencoder.bottleneck) <- rownames(filtered.counts)

amplicon.info <- amplicon.info[colnames(filtered.counts),]

methylome <- CreateSeuratObject(counts = t(filtered.counts), min.cells = 20, assay = "TAM" ) #only amplicons observed in at least 20 cells
methylome <- NormalizeData(methylome,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
methylome <- ScaleData(methylome, features = row.names(methylome))
methylome <- FindVariableFeatures(methylome,nfeatures = 200)
methylome <- RunPCA(methylome,npcs = 200)
methylome <- FindNeighbors(methylome, dims = 1:11)
methylome <- FindClusters(methylome, resolution = 0.05)
methylome <- RunUMAP(methylome, dims = 1:10)
DimPlot(methylome)
FeaturePlot(methylome, c('AMPL130986', 'AMPL130998', 'AMPL131004', 'AMPL134295'))
labels <- c("0"="Jurkat", "1" = "K562")
Idents(methylome) <- labels[as.character(Idents(methylome))]

methylome$log_count <- log10(methylome$nCount_TAM)

filtered.counts.uncut <- read.table("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/Sample2_80_percent/tsv/Sample2_80_percent.barcode.cell.distribution.tsv", row.names = 1, header=T)
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

allLikShort <- merge(allLikShort, amplicon.info[,c("amplicon", "Type")])
#allLikShort <- merge(allLikShort, bulk.methylation, all.x=T)
#compute dropout probability as a function of total read depth

#conclusion: False positive rate is near zero and fix, false negative rate can accurately be estimated
#fpr per amplicon
fpr <- subset(allLikShort, variable == "Jurkat" , select = c("amplicon", "uncut", "Type"))
rownames(fpr) <- fpr$amplicon
amplicon.info <- amplicon.info[rownames(methylome),]
fpr <- fpr[rownames(amplicon.info),]
#select relevant amplicons
filtered.counts <- filtered.counts[,rownames(methylome)]
amplicon.info$Type[is.na(amplicon.info$Type)] <- "NA"
selected <- filtered.counts[,(!grepl('Aci', fpr$Type)) & fpr$uncut > 0.9] > 0
selected <- apply(selected,2,as.numeric)
rownames(selected) <- rownames(filtered.counts)
non_hha_reads <- apply(filtered.counts[, row.names(amplicon.info[grepl('Aci', amplicon.info$Type),])], 1, sum)
rowinfo <- data.frame(row.names = rownames(filtered.counts), 
                      CellType = Idents(methylome), 
                      libsize = log10(methylome$nCount_TAM), 
                      NonCutAmplicons=log10(non_hha_reads), 
                      Nfeatures = methylome$nFeature_TAM,
                      doublet = ifelse(doublet_file[rownames(filtered.counts), 'DoubletDetectionLabel']==0, 'Singlet', 'Doublet'))

ph <- pheatmap(selected[row.names(rowinfo), ], 
               annotation_row = subset(rowinfo, select = c("CellType", "Nfeatures", "NonCutAmplicons")), 
               clustering_distance_cols = "binary", 
               clustering_distance_rows = "binary", 
               show_rownames = F, 
               show_colnames = F, 
               cutree_rows = 2, 
               clustering_method = "ward.D2",
               color=rev(inferno(50)),
               annotation_colors = color_map,
               fontsize=15)
ph$gtable$grobs[[1]]$gp <- gpar(lwd = 2)
ph$gtable$grobs[[2]]$gp <- gpar(lwd = 2)
png(file.path(plot_path, 'Supplemetary_CellLine_heatmap_selected.png'),
    width=1000,
    height=1600)
grid.newpage()
grid.draw(ph$gtable)
dev.off()
write.csv(rowinfo, '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Supplement/cell_metadata.csv')
write.csv(fpr[fpr$uncut > 0.9 & (!grepl('Aci', fpr$Type)),], '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Supplement/selected_amplicons.csv')
              