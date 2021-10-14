############## plot_denoised_UMAP.R ##############
#' Since the output of the autoencoder is a matrix of mixture proportions (values between 0 and 1), this matrix can also 
#' be visualized in two-dimensional space. Here, we show this for UMAPs and PCA plots. Additional inputs are
#' the original read count matrix, the cell assignment to cell types, and the amplicon information.

library(Seurat)
library(ggplot2)
library(pheatmap)
library(viridis)
sample <- 'BCells_Sample7_70_percent_good_performance'
meth.data <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/methylation_autoencoder/mixture_prob_dca.csv'),
                      header = TRUE,
                      row.names = 1)
input <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                    sep='\t',
                    header=TRUE)
mis.amplis <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/methylation_autoencoder/missing_amplicons.csv'),
                      header = TRUE,
                      row.names = 1)
input <- input[,!(colnames(input)%in%as.character(mis.amplis[,1]))]
ampli.info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.tsv')
cell.info <- read.csv('/users/mscherer/cluster/project/Methylome/infos/BCells/cell_sexes_Sample7.csv',row.names = 1)
colnames(meth.data) <- colnames(input)
row.names(meth.data) <- row.names(input)
seurat.obj <- CreateSeuratObject(t(meth.data),
                                 assay = "DNAm",
                                 meta.data = cell.info)
seurat.obj <- NormalizeData(seurat.obj,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
seurat.obj <- ScaleData(seurat.obj, features = row.names(seurat.obj))
seurat.obj <- FindVariableFeatures(seurat.obj,nfeatures = 700)
seurat.obj <- RunPCA(seurat.obj,npcs = 20)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:11)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:10)
DimPlot(seurat.obj, reduction = "umap")

DimPlot(seurat.obj, reduction = "umap", group.by = 'sex')

theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
FeaturePlot(seurat.obj,features = c('AMPL130835','AMPL130645'))
to.plot <- as.data.frame(seurat.obj[["umap"]]@cell.embeddings)
to.plot$AMPL130645 <- meth.data[row.names(to.plot),'AMPL130645'] #cg20821187
plot <- ggplot(to.plot,aes(x=UMAP_1,y=UMAP_2,color=AMPL130645))+geom_point()+theme
plot

ampli <- 'AMPL200226' #cg07303577
to.plot <- as.data.frame(seurat.obj[["umap"]]@cell.embeddings)
to.plot[,ampli] <- meth.data[row.names(to.plot),ampli]
plot <- ggplot(to.plot,aes_string(x='UMAP_1',y='UMAP_2',color=ampli))+geom_point()+theme
plot

ampli <- 'AMPL130620' #cg26653824
to.plot <- as.data.frame(seurat.obj[["umap"]]@cell.embeddings)
to.plot[,ampli] <- meth.data[row.names(to.plot),ampli]
plot <- ggplot(to.plot,aes_string(x='UMAP_1',y='UMAP_2',color=ampli))+geom_point()+theme
plot

to.plot <- as.data.frame(seurat.obj[["umap"]]@cell.embeddings)
to.plot <- data.frame(to.plo)