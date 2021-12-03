################### F2_cluster_UMAPs.R ################################
#' This file generates UMAPs of clusters and associated differential CpGs.

library(Seurat)
color_map <- list(CellType=c('naive'='#fcbd7e',
                             'memory1'='#fc6571',
                             'memory2'='#fc3262'),
                  Bulk=c('Naive B-cell high'='#fcbdbd',
                         'Memory B-cell high'='#fc00bd'),
                  Class=c('Naive B-cell high'='#fcbdbd',
                          'Memory B-cell high'='#fc00bd'),
                  CellType_reclustering=c('Cluster1'='#fcbd7e',
                                          'Cluster2a'='#fc3262',
                                          'Cluster2b'='#bd0071',
                                          'Cluster2c'='#8e008e'))
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=15),
                    axis.text=element_text(color='black',size=12),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    axis.text.y=element_blank())
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/cluster_umaps/'
rowinfo <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/rowinfo_reclustering_doublet.csv',
                    row.names = 1)
filtered.counts <- read.table("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/tsv/BCells_Sample7_70_percent_good_performance.barcode.cell.distribution_with_MCL.tsv", row.names = 1, header=T)
load('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/meth.data.numeric.Rdata')
cell_assignment <- read.table('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/MBC_assignment.txt')
meth.data.numeric <- meth.data.numeric[,as.character(cell_assignment$V5)]
for(clust1 in unique(rowinfo$CellType_reclustering)){
  for(clust2 in unique(rowinfo$CellType_reclustering)){
    diff_cpg_file <- paste0('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/differential/differential_CpGs_', clust1, 'vs', clust2, '.csv')
    if(!file.exists(diff_cpg_file)) next
    sel_data <- filtered.counts[row.names(rowinfo[rowinfo$CellType_reclustering%in%c(clust1, clust2), ]), diff_cpgs$Amplicon]
    seurat.obj <- CreateSeuratObject(t(sel_data),
                                     assay = "DNAm",
                                     meta.data = rowinfo)
    seurat.obj <- NormalizeData(seurat.obj,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000)
    seurat.obj <- ScaleData(seurat.obj, features = row.names(seurat.obj))
    seurat.obj <- FindVariableFeatures(seurat.obj,nfeatures = 200)
    seurat.obj <- RunPCA(seurat.obj,npcs = 200)
    ElbowPlot(seurat.obj)
    seurat.obj <- FindNeighbors(seurat.obj, dims = 1:11)
    seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:10)
    plot <- DimPlot(seurat.obj, reduction = 'umap', group.by = 'CellType_reclustering')+plot_theme+scale_color_manual(values=color_map[['CellType_reclustering']])
    ggsave(paste0(plot_path, 'UMAP_', clust1, 'vs', clust2, '.pdf'), plot)
  }
}
