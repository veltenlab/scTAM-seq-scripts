############## PCA_UMAP_Seurat.R ##############
#' This script creates UMAP and PCA plots for given read counts from scTAM-seq. Standard preprocessing as for
#' scRNA-seq data (log-normalizing, scaling, PCA, nearest neighbors) is executed. Further information about
#' cell clustering is visualized in the heatmap.

library(Seurat)
library(ggplot2)
sample <- 'Sample8_70_percent_good_performance'
dat <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                     header = T)
ampli_info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
rowinfo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/rowinfo.csv'),
                    row.names = 1)
dat <- dat[row.names(rowinfo), row.names(ampli_info)[ampli_info$Type.of.amplicon%in%'CpG.B.cell.diff']]
out.folder <- paste0('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/Sample8/')
cell.info <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/cell_sexes.csv'),
                      row.names = 1)
rowinfo$Donor <- NA 
rowinfo[row.names(cell.info), 'Donor'] <- cell.info[, 'sex']
seurat.obj <- CreateSeuratObject(t(dat),
                                 assay = "DNAm",
                                meta.data = rowinfo)
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
plot <- DimPlot(seurat.obj, reduction = "umap", group.by = 'CellType')+scale_color_manual(values=c('naive B-cells'='#fcbd7e',
                                                                                                     'memory B-cells'='#fc3262'))
ggsave('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/Sample8//UMAP.png', plot)
plot <- DimPlot(seurat.obj, reduction = "umap", group.by = 'Donor')#+scale_color_manual(values=c('naive B-cells'='#fcbd7e',
               #                                                                                    'memory B-cells'='#fc3262'))
ggsave('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/Sample8/UMAP_Donor.png', plot)

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

rowinfo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/rowinfo.csv'),
                    row.names = 1)
rowinfo <- subset(rowinfo, subset = CellType_broad=='memory B-cells')
dat <- dat[row.names(rowinfo), row.names(ampli_info)[ampli_info$Type.of.amplicon%in%'CpG.B.cell.diff']]
out.folder <- paste0('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/Sample8/')
cell.info <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/cell_sexes.csv'),
                      row.names = 1)
rowinfo$Donor <- NA 
rowinfo[row.names(cell.info), 'Donor'] <- cell.info[, 'sex']
seurat.obj <- CreateSeuratObject(t(dat),
                                 assay = "DNAm",
                                 meta.data = rowinfo)
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
plot <- DimPlot(seurat.obj, reduction = "umap", group.by = 'CellType')+scale_color_manual(values=c('cs-memory B-cells'='#8e008e',
                                                                                                   'ns-memory B-cells'='#fc3262'))
ggsave('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/Sample8//UMAP_cs_ns.png', plot)

FeaturePlot(seurat.obj, reduction = "umap", feature = 'CD10')
