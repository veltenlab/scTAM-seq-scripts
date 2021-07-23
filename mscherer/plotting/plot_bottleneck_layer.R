library(ggfortify)
sample <- 'Sample5_80_percent'
bottleneck.layer <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/methylation_autoencoder/bottleneck.csv'),
                             header=TRUE)
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
input <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                    sep='\t',
                    header=TRUE)
clust.file <- clust.file[row.names(input),]
clust.file <- clust.file[!(clust.file$CellType%in%'Mixed'),]
input <- input[row.names(clust.file),]
row.names(bottleneck.layer) <- row.names(input)
dat <- data.frame(bottleneck.layer,CellType=clust.file$CellType)
pca.obj <- prcomp(bottleneck.layer)
autoplot(pca.obj,data = dat, colour="CellType")

library(Seurat)
seurat.obj <- CreateSeuratObject(t(bottleneck.layer),
                                 assay = "DNAm",
                                 meta.data = clust.file)
seurat.obj <- NormalizeData(seurat.obj,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
seurat.obj <- ScaleData(seurat.obj, features = row.names(seurat.obj))
seurat.obj <- FindVariableFeatures(seurat.obj,nfeatures = 200)
seurat.obj <- RunPCA(seurat.obj,npcs = 5)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:5)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:5)
DimPlot(seurat.obj, reduction = "umap",group.by = 'CellType',
        cols = c('Jurkat'='#ff9289', 'K562'='#82b7ff','Mixed'='gray'))
