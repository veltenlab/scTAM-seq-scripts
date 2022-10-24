############## F2_C_UMAP.R ##############
#' This script creates UMAP and PCA plots for given read counts from scTAM-seq. Standard preprocessing as for
#' scRNA-seq data (log-normalizing, scaling, PCA, nearest neighbors) is executed. Further information about
#' cell clustering is visualized in the heatmap.

library(Signac)
library(Seurat)
library(ggplot2)
library(SeuratWrappers)
library(monocle3)
library(viridis)
library(gridExtra)
sample <- 'GSM5935918_Blood_HhaI'
cols <- c('naive B-cells'='#fcbd7e',
          '1'='#4a0f1d',
          '2_2'='#843c15',
          '2_1'='#c45a22',
          '4'='#f2615f',
          '5'='#fbe5d9',
          'cs-memory B-cells'='#8e008e')
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=7),
                    axis.text=element_text(color='black',size=6),
                    axis.title=element_blank(),
                    axis.ticks=element_line(color='black', size=.1),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.key=element_rect(color='black', fill=NA),
                    legend.position='none',
                    plot.title=element_blank())
plot_theme_facet <- theme(panel.background = element_rect(color='black',fill='white'),
                          panel.grid=element_blank(),
                          plot.background=element_blank(),
                          text=element_text(color='black',size=7),
                          axis.text=element_text(color='black',size=6),
                          axis.title=element_blank(),
                          axis.ticks=element_line(color='black', size=.1),
                          strip.background = element_blank(),
                          strip.text.x = element_blank(),
                          legend.position='none',
                          legend.key=element_rect(color='black', fill=NA),
                          legend.key.size = unit(.1, 'cm'), 
                          legend.key.height = unit(.1, 'cm'), 
                          legend.key.width = unit(.1, 'cm'), 
                          legend.title = element_text(size=7),
                          legend.text = element_text(size=6),
                          plot.title=element_text(color='black',size=6),
                          axis.ticks.length=unit(.1, "cm"))
dat <- read.table(paste0('../../data/',sample,'.tsv.gz'),
                  header = T)
rowinfo <- read.csv(paste0("../../misc/",sample,"/tsv/rowinfo_NCS.csv"),
                    row.names = 1)
rowinfo$ncBCellClustering_detailed <- rowinfo$ncBCellClustering
rowinfo$ncBCellClustering_detailed[rowinfo$ncBCellClustering_detailed%in%'Other'] <- rowinfo$CellType[rowinfo$ncBCellClustering_detailed%in%'Other']
colinfo <- read.csv(paste0("../../misc/GSM5935923_BM_undigested/tsv/colinfo.csv"),
                    row.names = 1)
out.folder <- '~'
dat <- dat[row.names(rowinfo), row.names(colinfo)]
dat <- ifelse(dat>0, 1, 0)
seurat.obj <- CreateSeuratObject(t(dat),
                                 assay = "DNAm",
                                 meta.data = rowinfo)
seurat.obj <- NormalizeData(seurat.obj,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
seurat.obj <- ScaleData(seurat.obj, features = row.names(seurat.obj))
seurat.obj <- RunPCA(seurat.obj,npcs = 200, features=row.names(seurat.obj))
ElbowPlot(seurat.obj)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:11)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:10)
DimPlot(seurat.obj, reduction = "umap")
plot <- DimPlot(seurat.obj, reduction = "umap")+plot_theme
ggsave(file.path(out.folder,"UMAP.png"), width = 40, height = 45, unit='mm')
to_plot <- as.data.frame(seurat.obj[['umap']]@cell.embeddings)
to_plot <- data.frame(to_plot, seurat.obj[[]])
to_plot <- to_plot[, c('UMAP_1',
                       'UMAP_2',
                       'ncBCellClustering')]
to_plot$ncBCellClustering[to_plot$ncBCellClustering%in%'Other'] <- seurat.obj[[]]$CellType[to_plot$ncBCellClustering%in%'Other']
plot <- ggplot(to_plot, aes(x=UMAP_1, y=UMAP_2, color=ncBCellClustering))+geom_point(size=.2, stroke=.2)+
  scale_color_manual(values=cols)+plot_theme
ggsave(file.path(out.folder,"UMAP_ncBCellClustering.png"), width = 40, height = 45, unit='mm')

to_plot <- as.data.frame(seurat.obj[['umap']]@cell.embeddings)
to_plot <- data.frame(to_plot, seurat.obj[[]])
to_plot <- to_plot[, c('UMAP_1',
                       'UMAP_2',
                       'CD11c',
                       'CD27')]
to_plot <- reshape2::melt(to_plot, id=c('UMAP_1', 'UMAP_2'))
to_plot$value <- ifelse(as.numeric(to_plot$value)>1, 'Expressed', 'NotExpressed')
plot_list <- list()
for(ab in unique(to_plot$variable)){
  sel_dat <- subset(to_plot, subset=variable==ab)
  plot <- ggplot(sel_dat, aes(x=UMAP_1, y=UMAP_2, color=value))+geom_point(size=.1, stroke=.1)+
    plot_theme_facet+ggtitle(ab)+scale_color_manual(values=c('Expressed'='#104e8b',
                                                             'NotExpressed'='gray80'))
  plot_list[[ab]] <- plot
}

png(file.path(out.folder,"UMAP_proteins.png"), width = 80, height = 28, units='mm', res=300)
plots <- do.call(grid.arrange, c(plot_list, ncol = 2))
dev.off()

Idents(seurat.obj) <- 'ncBCellClustering_detailed'
s1_cells <- names(Idents(seurat.obj)[Idents(seurat.obj)%in%'naive B-cells'])
seurat.mono <- as.cell_data_set(seurat.obj)
seurat.mono <- cluster_cells(cds = seurat.mono, reduction_method = "UMAP")
seurat.mono <- learn_graph(seurat.mono, use_partition = FALSE, close_loop=FALSE)

seurat.mono <- order_cells(seurat.mono, reduction_method = "UMAP", root_cells=s1_cells)

to_plot <- as.data.frame(seurat.obj[['umap']]@cell.embeddings)
to_plot <- data.frame(to_plot, Pseudotime=seurat.mono@principal_graph_aux$UMAP$pseudotime)
plot <- ggplot(to_plot, aes(x=UMAP_1, y=UMAP_2, color=Pseudotime))+geom_point(size=.2, stroke=.2)+
  plot_theme+scale_color_viridis(option='inferno')
ggsave(file.path(out.folder,"UMAP_pseudotime.png"), width = 40, height = 45, unit='mm')

