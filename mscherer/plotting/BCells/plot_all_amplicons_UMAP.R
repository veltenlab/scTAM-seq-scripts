############## plot_all_amplicons_UMAP.R ##############

#.libPaths(c(.libPaths(),'/users/mscherer/conda/envs/rnbeads/lib/R/library/'))
library(Seurat)
library(ggplot2)
library(pheatmap)
library(viridis)
library(RnBeads)
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
colnames(meth.data) <- colnames(input)
row.names(meth.data) <- row.names(input)
seurat.obj <- CreateSeuratObject(t(meth.data),
                                 assay = "DNAm",
                                 meta.data = input)
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
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
cpg.names <- row.names(rnb.annotation2data.frame(rnb.get.annotation('probesEPIC')))
anno <- unlist(rnb.get.annotation('probesEPIC'))               
ampli.info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.tsv')
plot.path <- paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/plots/amplicons_wo_MCL/')
plot.raw <- paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/plots/amplicons_wo_MCL/raw/')
for(ampli in colnames(meth.data)){
  chr.code <- paste0(ampli.info[ampli,'chr'],
                     ':',
                     ampli.info[ampli,'amplicon_start'],
                     '-',
                     ampli.info[ampli,'amplicon_end'])
  my.gr <- GRanges(chr.code)
  op <- findOverlaps(my.gr,anno)
  cpg <- paste(cpg.names[subjectHits(op)],collapse='_')
  to.plot <- as.data.frame(seurat.obj[["umap"]]@cell.embeddings)
  to.plot[,ampli] <- meth.data[row.names(to.plot),ampli]
  plot <- ggplot(to.plot,aes_string(x='UMAP_1',y='UMAP_2',color=ampli))+geom_point()+theme+
    scale_color_gradientn(colors=rev(viridis(15)),
                          limits=c(0,1))
  f.name <- file.path(plot.path,paste0(ampli,'_',cpg,'.pdf'))
  ggsave(f.name,plot)
  to.plot <- as.data.frame(seurat.obj[["umap"]]@cell.embeddings)
  to.plot[,ampli] <- log10(input[row.names(to.plot),ampli]+1)
  plot <- ggplot(to.plot,aes_string(x='UMAP_1',y='UMAP_2',color=ampli))+geom_point()+theme+
    scale_color_gradientn(colors=rev(magma(15)), name='log10(reads)')
  f.name <- file.path(plot.raw,paste0(ampli,'_',cpg,'.pdf'))
  ggsave(f.name,plot)
}
