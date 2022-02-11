################### Supplement_comparison_with_Triana_et_al.R ################### 
#' With this file, we compare the results we obtained (expression of different surface markers)
#' with the results obtained in Triana et al, Nat Imm, 2021

library(Seurat)
seurat.obj <- readRDS('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Sergio/data/Healthy.rds')
# only the aged samples or also include the BM cells?
seurat.obj <- subset(seurat.obj, idents='Nonswitched memory B cells')
seurat.obj <- FindVariableFeatures(seurat.obj)
seurat.obj <- RunPCA(seurat.obj, reduction.name = 'PCA_memory')
seurat.obj <- FindNeighbors(seurat.obj, reduction = 'PCA_memory', graph.name = 'KNN_memory')
seurat.obj <- FindClusters(seurat.obj, graph.name = 'KNN_memory', resolution = 0.1)
seurat.obj <- RunUMAP(seurat.obj, reduction='PCA_memory', reduction.name='UMAP_memory', dims=1:10)
DimPlot(seurat.obj, reduction = 'UMAP_memory')
FeaturePlot(seurat.obj, c('CD27-AB', 'CD1c-AB'), reduction = 'UMAP_memory')
FeaturePlot(seurat.obj, 'CD45-AB', reduction = 'UMAP_memory')
Idents(seurat.obj) <- 'KNN_memory_res.0.1'
diff.markers <- FindMarkers(seurat.obj, ident.1 = '0', ident.2 = '1')
head(diff.markers, 10) 

diff.markers <- FindMarkers(seurat.obj, ident.1 = '1', ident.2 = '2')
head(diff.markers, 10) 

library(Seurat)
library(ggplot2)
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_line(color='black', size=.25),
                    strip.background = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5),
                    legend.position='none')
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Supplement/'
seurat.obj <- readRDS('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Sergio/data/Healthy.rds')
seurat.obj <-seurat.obj[, Idents(seurat.obj)=='Nonswitched memory B cells'|Idents(seurat.obj)=='Class switched memory B cells']
to_plot <- GetAssayData(seurat.obj)[c('CD27-AB',
                                      'CD25-AB',
                                      'CD1c-AB'), ]
to_plot <- data.frame(t(to_plot), CellType=Idents(seurat.obj))
to_plot <- reshape2::melt(to_plot, id='CellType')
colnames(to_plot) <- c('CellType', 'Antibody', 'Expression')
plot <- ggplot(to_plot, aes(x=CellType, y=Expression, fill=CellType))+geom_violin()+geom_boxplot(alpha=0.25)+
  facet_wrap(Antibody~., ncol=2)+
  plot_theme+scale_fill_manual(values=c('Nonswitched memory B cells'='#fc3262',
                                        'Class switched memory B cells'='#8e008e'))
ggsave(file.path(plot_path, 'Triana_memory_boxplots.pdf'), plot, width=200, height=75, unit='mm')
