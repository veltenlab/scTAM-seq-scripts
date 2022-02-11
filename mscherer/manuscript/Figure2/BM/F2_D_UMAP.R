############## F2_D_UMAP.R ##############
#' This script creates UMAP and PCA plots for given read counts from scTAM-seq. Standard preprocessing as for
#' scRNA-seq data (log-normalizing, scaling, PCA, nearest neighbors) is executed. Further information about
#' cell clustering is visualized in the heatmap.

library(Signac)
library(Seurat)
library(ggplot2)
library(SeuratWrappers)
library(monocle3)
sample <- 'Sample11_70_percent_good_performance'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=15),
                    axis.text=element_text(color='black',size=15),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.key=element_rect(color='black', fill=NA))
dat <- read.table(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                     header = T)
rowinfo <- read.csv(paste0("/users/lvelten/project/Methylome/analysis/missionbio/BM/",sample,"/tsv/rowinfo.csv"),
                    row.names = 1)
#colinfo <- read.table('/users/lvelten/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt',
#                      row.names = 8)
colinfo <- read.csv(paste0("/users/lvelten/project/Methylome/analysis/missionbio/BM/",sample,"/tsv/colinfo.csv"),
                    row.names = 1)
out.folder <- '/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure2/Sample11/'
dat <- dat[row.names(rowinfo), row.names(colinfo)]
#dat <- dat[row.names(rowinfo), ]
dat <- ifelse(dat>0, 1, 0)
seurat.obj <- CreateSeuratObject(t(dat),
                                 assay = "DNAm",
                                meta.data = rowinfo)
seurat.obj <- NormalizeData(seurat.obj,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
seurat.obj <- ScaleData(seurat.obj, features = row.names(seurat.obj))
#seurat.obj <- FindVariableFeatures(seurat.obj,nfeatures = 200)
seurat.obj <- RunPCA(seurat.obj,npcs = 200, features=row.names(seurat.obj))
ElbowPlot(seurat.obj)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:11)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:10)
DimPlot(seurat.obj, reduction = "umap")
plot <- DimPlot(seurat.obj, reduction = "umap")+plot_theme
ggsave(file.path(out.folder,"UMAP.png"), width = 125, height = 100, unit='mm')
plot <- DimPlot(seurat.obj, reduction = "umap", group.by='CellType_detailed')+plot_theme+scale_color_manual(values=c('naive B-cells'='#fcbd7e',
                                                                                                                     'memory B-cells'='#fc6571',
                                                                                                                     'S2 cells'='#d7ef7e',
                                                                                                                     'S3-S4 cells'='#fcef7e',
                                                                                                                      'S1 cells'='#bdfc7e'))
ggsave(file.path(out.folder,"UMAP_CellType.png"), width = 125, height = 100, unit='mm')

#seurat.obj <- RunTSNE(seurat.obj, dims = 1:10)
#plot <- DimPlot(seurat.obj, reduction = "tsne", group.by='CellType')+plot_theme+scale_color_manual(values=c('naive B-cells'='#fcbd7e',
#                                                                                                            'memory B-cells'='#fc6571',
#                                                                                                            'S2-S4 cells'='#fcef7e',
#                                                                                                            'S1 cells'='#bdfc7e'))
#ggsave(file.path(out.folder,"tSNE_CellType.png"), width = 125, height = 100, unit='mm')



#diff.markers <- FindAllMarkers(seurat.obj)

Idents(seurat.obj) <- 'CellType'
s1_cells <- names(Idents(seurat.obj)[Idents(seurat.obj)%in%'S1 cells'])
seurat.mono <- as.cell_data_set(seurat.obj)
seurat.mono <- cluster_cells(cds = seurat.mono, reduction_method = "UMAP")
seurat.mono <- learn_graph(seurat.mono, use_partition = FALSE, close_loop=FALSE)

seurat.mono <- order_cells(seurat.mono, reduction_method = "UMAP", root_cells=s1_cells)

# plot trajectories colored by pseudotime
plot <- plot_cells(
  cds = seurat.mono,
  color_cells_by = "pseudotime",
  label_branch_points = FALSE,
  label_leaves=FALSE,
  label_roots=FALSE
)+plot_theme
ggsave(file.path(out.folder,"UMAP_trajectory.png"), width = 125, height = 100, unit='mm')
