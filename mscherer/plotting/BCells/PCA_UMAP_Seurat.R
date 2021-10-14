############## PCA_UMAP_Seurat.R ##############
#' This script creates UMAP and PCA plots for given read counts from scTAM-seq. Standard preprocessing as for
#' scRNA-seq data (log-normalizing, scaling, PCA, nearest neighbors) is executed. Further information about
#' cell clustering is visualized in the heatmap.

library(Seurat)
library(ggplot2)
sample <- 'BCells_Sample7_70_percent_good_performance'
dat <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                     header = T)
out.folder <- paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/plots/')
seurat.obj <- CreateSeuratObject(t(dat),
                                 assay = "DNAm")
seurat.obj <- NormalizeData(seurat.obj,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
seurat.obj <- ScaleData(seurat.obj, features = row.names(seurat.obj))
seurat.obj <- FindVariableFeatures(seurat.obj,nfeatures = 700)
seurat.obj <- RunPCA(seurat.obj,npcs = 20)
ElbowPlot(seurat.obj)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:11)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:10)
DimPlot(seurat.obj, reduction = "umap")
plot <- DimPlot(seurat.obj, reduction = "umap")
ggsave(file.path(out.folder,"UMAP.pdf"),width = 8, height = 5)
plot <- DimPlot(seurat.obj, reduction = "pca")
ggsave(file.path(out.folder,"PCA.pdf"),width = 8, height = 5)

FeaturePlot(seurat.obj,features = c('AMPL121673','AMPL130599'))
diff.res <- FindMarkers(seurat.obj, ident.1 = 1)
FeaturePlot(seurat.obj,features = c('AMPL130835','AMPL130645'))

theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
cut.off <- 5
to.plot <- as.data.frame(seurat.obj[["umap"]]@cell.embeddings)
to.plot$AMPL130645 <- as.numeric(dat[row.names(to.plot),'AMPL130645']>cut.off)
plot <- ggplot(to.plot,aes(x=UMAP_1,y=UMAP_2,color=AMPL130645))+geom_point()+theme
plot
