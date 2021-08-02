############## PCA_UMAP_Seurat.R ##############
#' This script creates UMAP and PCA plots for given read counts from scTAM-seq. Standard preprocessing as for
#' scRNA-seq data (log-normalizing, scaling, PCA, nearest neighbors) is executed. Further information about
#' cell clustering is visualized in the heatmap.

library(Seurat)
library(ggplot2)
sample <- 'Sample5_80_percent'
dat <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                     header = T)
clust.file <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/tsv/cluster_assignment.tsv"))
out.folder <- paste0('/users/mscherer/cluster/project/Methylome/analysis/Seurat/',sample)
dat <- dat[row.names(clust.file),]
dat <- dat[clust.file$CellType!="Mixed",]
clust.file <- clust.file[clust.file$CellType!="Mixed",]
dat <- dat[clust.file$Barcode,]
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
                cols = c('Jurkat'='#ff9289', 'K562'='#82b7ff','Mixed'='gray'))
ggsave(file.path(out.folder,"UMAP.pdf"),width = 8, height = 5)
plot <- DimPlot(seurat.obj, reduction = "pca",group.by = 'CellType',
                cols = c('Jurkat'='#ff9289', 'K562'='#82b7ff','Mixed'='gray'))
ggsave(file.path(out.folder,"PCA.pdf"),width = 8, height = 5)
