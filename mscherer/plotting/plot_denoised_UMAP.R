library(Seurat)
library(ggplot2)
library(pheatmap)
library(viridis)
sel.amplis <- FALSE
drop.cells <- FALSE
sample <- 'Sample2_80_percent'
meth.data <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/methylation_autoencoder/mixture_prob_dca.csv'),
                      header = TRUE,
                      row.names = 1)
input <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                    sep='\t',
                    header=TRUE)

is.zero.inpu <- apply(input,2,function(x)all(x==0))
input <- input[,!is.zero.inpu]
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
if(drop.cells){
  rem.cells <- clust.file$Barcode[clust.file$CellType=='Mixed']
  input <- input[-which(row.names(input)%in%rem.cells),]
}
ampli.info <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/amplicons/cell_lines/R/amplicon_info_with_RnBeads.csv')
ampli.info <- ampli.info[!is.zero.inpu,]
row.names(ampli.info) <- ampli.info$AmpID_design_1898
if(sel.amplis){
  rel.info <- read.csv('/users/mscherer/cluster/project/Methylome/infos/cell_lines/analysis_Agostina_Apr21/Summary_methylation_values_amplicons_K562_Jurkat.csv')
  high.conf <- rel.info$AmpID_design_1898[grepl('Jurkat.low.K562.high|Jurkat.high.K562.low',rel.info$Label)]
  ampli.info <- ampli.info[high.conf,]
  input <- input[,high.conf]
}
colnames(meth.data) <- colnames(input)
row.names(meth.data) <- row.names(input)
clust.file <- clust.file[row.names(meth.data),]
row.names(clust.file) <- row.names(meth.data)
seurat.obj <- CreateSeuratObject(t(meth.data),
                                 assay = "DNAm",
                                meta.data = clust.file)
seurat.obj <- NormalizeData(seurat.obj,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
seurat.obj <- ScaleData(seurat.obj, features = row.names(seurat.obj))
seurat.obj <- FindVariableFeatures(seurat.obj,nfeatures = 200)
seurat.obj <- RunPCA(seurat.obj,npcs = 200)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:11)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:10)
DimPlot(seurat.obj, reduction = "umap",group.by = 'CellType',
                cols = c('Jurkat'='#ff9289', 'K562'='#82b7ff','Mixed'='gray'))
