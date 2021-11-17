############## Supplement_Cell_Lines_PCA_UMAP.R ##############
#' This script creates UMAP and PCA plots for given read counts from scTAM-seq. Standard preprocessing as for
#' scRNA-seq data (log-normalizing, scaling, PCA, nearest neighbors) is executed. Further information about
#' cell clustering is visualized in the heatmap.

library(Seurat)
library(ggplot2)
sample <- 'Sample5_80_percent'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=15),
                    axis.text=element_text(color='black',size=15),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.key=element_rect(color='black', fill=NA))
dat <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                     header = T)
clust.file <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/tsv/cluster_assignment.tsv"))
out.folder <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Supplement/'
dat <- dat[row.names(clust.file),]
dat <- dat[clust.file$Barcode,]
doublet.file <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores_DoubletDetection.csv'),
                         row.names=2)
clust.file$DoubletDetectionScore <- doublet.file[row.names(clust.file), 'DoubletDetectionScore']
seurat.obj <- CreateSeuratObject(t(dat),
                                 assay = "DNAm",
                                meta.data = clust.file)
seurat.obj <- NormalizeData(seurat.obj,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
seurat.obj <- ScaleData(seurat.obj, features = row.names(seurat.obj))
seurat.obj <- FindVariableFeatures(seurat.obj,nfeatures = 200)
seurat.obj <- RunPCA(seurat.obj,npcs = 200)
ElbowPlot(seurat.obj)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:11)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:10)
DimPlot(seurat.obj, reduction = "umap")
plot <- DimPlot(seurat.obj, reduction = "umap",group.by = 'CellType',
                cols = c('Jurkat'='#ff9289', 'K562'='#82b7ff','Mixed'='gray'))+plot_theme
ggsave(file.path(out.folder,"Supplement_UMAP_with_doublets.pdf"), width = 125, height = 100, unit='mm')
plot <- FeaturePlot(seurat.obj, reduction = "umap",features = 'DoubletDetectionScore')+plot_theme
ggsave(file.path(out.folder,"Supplement_UMAP_DoubletDetectionScore.pdf"), width = 125, height = 100, unit='mm')
