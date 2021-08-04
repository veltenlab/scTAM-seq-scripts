library(Seurat)
library(ggplot2)
library(DoubletFinder)
sample <- 'Sample5_80_percent'
dat <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                  header = T)
clust.file <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/tsv/cluster_assignment.tsv"))
out.folder <- paste0('/users/mscherer/cluster/project/Methylome/analysis/Seurat/',sample)
dat <- dat[row.names(clust.file),]
dat <- dat[clust.file$CellType!="Mixed",]
clust.file <- clust.file[clust.file$CellType!="Mixed",]
dat <- dat[clust.file$Barcode,]
seurat.obj <- CreateSeuratObject(raw.data=t(dat),
                                  counts=t(dat),
                                 assay = "DNAm",
                                 meta.data = clust.file)
seurat.obj <- NormalizeData(seurat.obj,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
seurat.obj <- ScaleData(seurat.obj, features = row.names(seurat.obj))
seurat.obj <- FindVariableFeatures(seurat.obj,nfeatures = 200)
seurat.obj <- RunPCA(seurat.obj,npcs = 200)
seurat.obj@assays$DNAm@raw.data <- dat
seurat.obj <- doubletFinder(seurat.obj)
