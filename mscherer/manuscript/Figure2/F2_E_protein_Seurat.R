############## PCA_UMAP_Seurat.R ##############
#' This script creates UMAP and PCA plots for given read counts from scTAM-seq. Standard preprocessing as for
#' scRNA-seq data (log-normalizing, scaling, PCA, nearest neighbors) is executed. Further information about
#' cell clustering is visualized in the heatmap.

library(Seurat)
library(Matrix)
library(ggplot2)
color_map <- list(CellType_broad=c('naive B-cells'='#fcbd7e',
                                   'memory B-cells'='#fc6571',
                                   'memory B-cells'='#fc3262'),
                  Bulk=c('Naive B-cell high'='#fcbdbd',
                         'Memory B-cell high'='#fc00bd'),
                  Class=c('Naive B-cell high'='#fcbdbd',
                          'Memory B-cell high'='#fc00bd'),
                  CellType=c('naive B-cells'='#fcbd7e',
                             'ns-memory B-cells'='#fc3262',
                             'cs-memory B-cells'='#8e008e'))
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
sample <- 'Sample11_70_percent_good_performance'
protein <- 'Sample11'
dat <- read.table(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                     header = T)
ampli_info <- read.table('/users/lvelten/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
rowinfo <- read.csv(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/', sample, '/tsv/rowinfo.csv'),
                    row.names = 1)
#rowinfo <- rowinfo[, 1:10]
dat <- dat[row.names(rowinfo), row.names(ampli_info)[ampli_info$Type.of.amplicon%in%'CpG.B.cell.diff']]
out.folder <- paste0('/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure2/Sample11/')
prot_data_all <- read.table(paste0('/users/lvelten/project/Methylome/analysis/missionbio/tapestri/protein/', 
                               protein, '/pipeline/results/tsv/Sample11-counts.tsv'),
                        header=TRUE)
prot_data_all <- prot_data_all[prot_data_all$cell_barcode%in%row.names(rowinfo), ]
prot_data <- matrix(nrow = nrow(rowinfo), ncol=length(unique(prot_data_all$ab_description)))
row.names(prot_data) <- row.names(rowinfo)
colnames(prot_data) <- unique(prot_data_all$ab_description)
for(ab in unique(prot_data_all$ab_description)){
  sel_rows <- subset(prot_data_all, subset=ab_description==ab)
  prot_data[sel_rows$cell_barcode, ab] <- sel_rows$raw
}
prot_data[is.na(prot_data)] <- 0
seurat.obj <- CreateSeuratObject(t(prot_data),
                                 assay = "SurfaceAntibodies",
                                meta.data = rowinfo)
seurat.obj <- NormalizeData(seurat.obj,
                            normalization.method = "CLR")
seurat.obj <- ScaleData(seurat.obj, features = row.names(seurat.obj))
#seurat.obj <- FindVariableFeatures(seurat.obj, nfeatures = 25)
seurat.obj <- RunPCA(seurat.obj, features=row.names(seurat.obj))
# seurat.obj <- RunPCA(seurat.obj, features=c('CD27',
#                                             'CD1c',
#                                             'CD71',
#                                             'CD69',
#                                             'CD62L',
#                                             'CD141',
#                                             'CD11c',
#                                             'CD22',
#                                             'HLA.DR',
#                                             'CD33',
#                                             'CD49d',
#                                             'CD5',
#                                             'CD11b',
#                                             'CD25',
#                                             'CD13'))
ElbowPlot(seurat.obj)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:5)
#seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)
seurat.obj <- RunUMAP(seurat.obj, dims=1:5)
plot <- DimPlot(seurat.obj, reduction = "umap", group.by = 'CellType_detailed')+scale_color_manual(values=c('naive B-cells'='#fcbd7e',
                                                                                                   'memory B-cells'='#fc6571',
                                                                                                   'S2 cells'='#d7ef7e',
                                                                                                   'S3-S4 cells'='#fcef7e',
                                                                                                   'S1 cells'='#bdfc7e'))
ggsave('/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure2/Sample11/protein_UMAP.png', plot)
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
ggsave('/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure2/Sample11/protein_UMAP.png', plot, width = 40, height = 45, unit='mm')

plot <- FeaturePlot(seurat.obj, reduction = "umap", feature = 'CD27')
ggsave('/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure2/Sample11/protein_UMAP_CD27.png', plot)

for_rowinfo <- t(GetAssayData(seurat.obj))
for_rowinfo_raw <- t(GetAssayData(seurat.obj, slot='counts'))
colnames(for_rowinfo_raw) <- paste0(colnames(for_rowinfo), '_raw')
rowinfo <- data.frame(rowinfo, for_rowinfo, for_rowinfo_raw)
write.csv(rowinfo, paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/', sample, '/tsv/rowinfo.csv'))
