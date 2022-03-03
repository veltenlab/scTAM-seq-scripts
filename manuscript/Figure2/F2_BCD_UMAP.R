############## F2_B_D_UMAP.R ##############
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
sample <- 'Sample11_70_percent_good_performance'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
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
                          text=element_text(color='black',size=6),
                          axis.text=element_text(color='black',size=5),
                          axis.title=element_blank(),
                          axis.ticks=element_line(color='black', size=.1),
                          strip.background = element_blank(),
                          strip.text.x = element_blank(),
                          legend.position='none',
                          legend.key=element_rect(color='black', fill=NA),
                          legend.key.size = unit(.1, 'cm'), 
                          legend.key.height = unit(.1, 'cm'), 
                          legend.key.width = unit(.1, 'cm'), 
                          legend.title = element_text(size=6),
                          legend.text = element_text(size=5),
                          plot.title=element_text(color='black',size=5),
                          axis.ticks.length=unit(.1, "cm"))
dat <- read.table(paste0('../../data/', sample, '.barcode.cell.distribution.tsv'),
                  header = T)
rowinfo <- read.csv(paste0("../../misc/", sample, "/tsv/rowinfo.csv"),
                    row.names = 1)
colinfo <- read.csv(paste0("../../misc/",sample,"/tsv/colinfo.csv"),
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
                       'CellType_detailed')]
plot <- ggplot(to_plot, aes(x=UMAP_1, y=UMAP_2, color=CellType_detailed))+geom_point(size=.2, stroke=.2)+
  plot_theme+scale_color_manual(values=c('naive B-cells'='#fcbd7e',
                                         'memory B-cells'='#fc6571',
                                         'S2 cells'='#d7ef7e',
                                         'S3-S4 cells'='#fcef7e',
                                         'S1 cells'='#bdfc7e'))
ggsave(file.path(out.folder,"UMAP_CellType.png"), width = 40, height = 45, unit='mm')


to_plot <- as.data.frame(seurat.obj[['umap']]@cell.embeddings)
to_plot <- data.frame(to_plot, seurat.obj[[]])
to_plot <- to_plot[, c('UMAP_1',
                       'UMAP_2',
                       'CD10',
                       'CD19',
                       'CD27',
                       'CD34',
                       'CD38',
                       'CD45RA',
                       'CD90'                       )]
to_plot <- reshape2::melt(to_plot, id=c('UMAP_1', 'UMAP_2'))
to_plot$value <- as.numeric(to_plot$value)
plot_list <- list()
for(ab in unique(to_plot$variable)){
  sel_dat <- subset(to_plot, subset=variable==ab)
  plot <- ggplot(sel_dat, aes(x=UMAP_1, y=UMAP_2, color=value))+geom_point(size=.1, stroke=.1)+
    plot_theme_facet+scale_color_viridis(option='mako', direction=-1)+ggtitle(ab)
  plot_list[[ab]] <- plot
}
png(file.path(out.folder,"UMAP_proteins.png"), width = 170, height = 40, units='mm', res=300)
plots <- do.call(grid.arrange, c(plot_list, ncol = 7))
dev.off()

to_plot <- as.data.frame(seurat.obj[['umap']]@cell.embeddings)
to_plot <- data.frame(to_plot, seurat.obj[[]])
to_plot <- to_plot[, c('UMAP_1',
                       'UMAP_2',
                       'CD138',
                       'CD5')]
to_plot <- reshape2::melt(to_plot, id=c('UMAP_1', 'UMAP_2'))
to_plot$value <- as.numeric(to_plot$value)
plot_list <- list()
for(ab in unique(to_plot$variable)){
  sel_dat <- subset(to_plot, subset=variable==ab)
  plot <- ggplot(sel_dat, aes(x=UMAP_1, y=UMAP_2, color=value))+geom_point(size=.15, stroke=.15)+
    plot_theme_facet+scale_color_viridis(option='mako', direction=-1)+ggtitle(ab)
  plot_list[[ab]] <- plot
}
png(file.path(out.folder,"UMAP_proteins_CD138_CD5.png"), width = 50, height = 50, units='mm', res=300)
plots <- do.call(grid.arrange, c(plot_list))
dev.off()

Idents(seurat.obj) <- 'CellType'
s1_cells <- names(Idents(seurat.obj)[Idents(seurat.obj)%in%'S1 cells'])
seurat.mono <- as.cell_data_set(seurat.obj)
seurat.mono <- cluster_cells(cds = seurat.mono, reduction_method = "UMAP")
seurat.mono <- learn_graph(seurat.mono, use_partition = FALSE, close_loop=FALSE)

seurat.mono <- order_cells(seurat.mono, reduction_method = "UMAP", root_cells=s1_cells)

# plot trajectories colored by pseudotime
to_plot <- as.data.frame(seurat.obj[['umap']]@cell.embeddings)
to_plot <- data.frame(to_plot, Pseudotime=seurat.mono@principal_graph_aux$UMAP$pseudotime)
plot <- ggplot(to_plot, aes(x=UMAP_1, y=UMAP_2, color=Pseudotime))+geom_point(size=.2, stroke=.2)+
  plot_theme+scale_color_viridis(option='inferno')
ggsave(file.path(out.folder,"UMAP_trajectory.png"), width = 40, height = 45, unit='mm')
