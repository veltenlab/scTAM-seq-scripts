############## F2_integration.R ##############
#' This script employs integration between scRNA-seq and scTAM-seq data

library(Signac)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(SeuratWrappers)
library(monocle3)
library(viridis)
library(gridExtra)
library(RnBeads)
sample <- 'Sample11_70_percent_good_performance'
protein <- 'Sample11'
out.folder <- '/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure2/Sample11/'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.title=element_blank(),
                    axis.ticks=element_line(color='black', size=.1),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    #legend.key=element_rect(color='black', fill=NA),
                    legend.key.size = unit(2, 'mm'), #change legend key size
                    legend.key.height = unit(2, 'mm'), #change legend key height
                    legend.key.width = unit(2, 'mm'), #change legend key width
                    legend.title = element_text(color='black',size=6), #change legend title font size
                    legend.text = element_text(color='black',size=5),
                    legend.position='none',
                    plot.title=element_blank())
dat <- read.table(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                  header = T)
rowinfo <- read.csv(paste0("/users/lvelten/project/Methylome/analysis/missionbio/BM/",sample,"/tsv/rowinfo.csv"),
                    row.names = 1)
colinfo <- read.csv(paste0("/users/lvelten/project/Methylome/analysis/missionbio/BM/",sample,"/tsv/colinfo.csv"),
                    row.names = 1)
out.folder <- '/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Figure2/Sample11/'
dat <- dat[row.names(rowinfo), row.names(colinfo)]
dat <- ifelse(dat>0, 1, 0)
seurat.obj <- CreateSeuratObject(t(dat),
                                 assay = "DNAm",
                                 meta.data = rowinfo)
dat <- read.table(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                     header = T)
ampli_info <- read.table('/users/lvelten/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
rowinfo <- read.csv(paste0('/users/lvelten/project/Methylome/analysis/missionbio/BM/', sample, '/tsv/rowinfo.csv'),
                    row.names = 1)
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
colnames(prot_data) <- paste0(colnames(prot_data), '-AB')
seurat.obj <- NormalizeData(seurat.obj,
                            normalization.method = "LogNormalize")
seurat.obj <- ScaleData(seurat.obj, features = row.names(seurat.obj))
adt_assay <- CreateAssayObject(counts = t(prot_data))
seurat.obj[["AB"]] <- adt_assay
DefaultAssay(seurat.obj) <- 'AB'
seurat.obj <- NormalizeData(seurat.obj,
                            normalization.method = "CLR")
seurat.obj <- ScaleData(seurat.obj,
    features = row.names(seurat.obj))
seurat.obj <- RunPCA(seurat.obj,
    features=row.names(seurat.obj))

sergio.obj <- readRDS('/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Sergio/data/Healthy.rds')
cells <- Idents(sergio.obj)
DefaultAssay(sergio.obj) <- 'RNA'
sergio.obj[['BOTH']] <- NULL
sergio.obj[['integrated']] <- NULL
sergio.obj <- FindVariableFeatures(sergio.obj)
sergio.obj <- RunPCA(sergio.obj,
 reduction.name='pcaRNA')
DefaultAssay(sergio.obj) <- 'AB'
VariableFeatures(sergio.obj) <- rownames(sergio.obj[["AB"]])
sergio.obj <- RunPCA(sergio.obj,
 reduction.name='pcaAB')
sergio.obj <- sergio.obj[, which(grepl('B cell', cells)|grepl('HSC', cells))]
sergio.obj <- FindMultiModalNeighbors(
  sergio.obj,
  reduction.list = list("pcaAB", "pcaRNA"), 
  dims.list = list(1:10, 1:10),
  modality.weight.name = "RNA.weight"
)
sergio.obj <- RunSPCA(sergio.obj,
    assay = 'AB',
    graph = 'wsnn',
    npcs = 10)
joint_features <- intersect(rownames(sergio.obj), rownames(seurat.obj))
anchors <- FindTransferAnchors(
  reference = sergio.obj,
  query = seurat.obj,
  normalization.method = "LogNormalize",
  reference.reduction = "spca",
  dims = 1:10,
  features = joint_features
)
sergio.obj$TrianaCellType <- Idents(sergio.obj)
seurat.obj <- MapQuery(
  anchorset = anchors,
  query = seurat.obj,
  reference = sergio.obj,
  refdata = list(
    ct = "TrianaCellType"
  ),
  reference.reduction = "spca"
)

DefaultAssay(seurat.obj) <- 'DNAm'
seurat.obj <- RunPCA(seurat.obj, npcs = 200, features=row.names(seurat.obj)[!grepl('-AB', row.names(seurat.obj))])
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:11)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:10)
to_plot <- as.data.frame(seurat.obj[['umap']]@cell.embeddings)
to_plot <- data.frame(to_plot, seurat.obj[[]])
to_plot <- to_plot[, c('UMAP_1',
                       'UMAP_2',
                       'predicted.ct')]
color_map <- c('Mature naive B cells'='#fcbd7e',
    'Nonswitched memory B cells'='#fc3262',
    'Class switched memory B cells'='#8e008e',
    'Immature B cells'='#fcef7e',
    'Pre-pro-B cells'='#d7ef7e',
    'Pro-B cells'='#d7ef7e',
    'Pre-B cells'='#d7ef7e',
    'HSCs & MPPs'='#bdfc7e')
plot <- ggplot(to_plot, aes(x=UMAP_1, y=UMAP_2, color=predicted.ct))+geom_point(size=.2, stroke=.2)+
  plot_theme+scale_color_manual(values=color_map)+ guides(color = guide_legend(override.aes = list(size = 1.5)))
ggsave(file.path(out.folder,"UMAP_Triana_reference.png"), width = 40, height = 45, unit='mm')

seurat.obj$CT_transformed <- c('S1 cells'='HSCs & MPPs',
    'S2 cells'='Pre-B cells',
    'S3-S4 cells'='Immature B cells',
    'memory B-cells'='Nonswitched memory B cells',
    'naive B-cells'='Mature naive B cells'
)[seurat.obj$CellType_detailed]
Idents(seurat.obj) <- 'CT_transformed'
sergio.obj$id <- 'reference'
seurat.obj$id <- 'query'
joint.obj <- merge(sergio.obj, seurat.obj)
joint.obj[["spca"]] <- merge(sergio.obj[["spca"]], seurat.obj[["ref.spca"]])
joint.obj <- RunUMAP(joint.obj,
    reduction = 'spca',
    dims = 1:10)
to_plot <- as.data.frame(joint.obj[['umap']]@cell.embeddings)
to_plot <- data.frame(to_plot, Idents=Idents(joint.obj))
to_plot <- to_plot[, c('UMAP_1',
                       'UMAP_2',
                       'Idents')]
plot <- ggplot(to_plot, aes(x=UMAP_1, y=UMAP_2, color=Idents))+geom_point(size=.2, stroke=.2)+
  plot_theme+scale_color_manual(values=color_map)+ guides(color = guide_legend(override.aes = list(size = 1.5)))
ggsave(file.path(out.folder,"UMAP_protein_labels.png"), width = 40, height = 45, unit='mm')

DefaultAssay(seurat.obj) <- 'DNAm'
DefaultAssay(joint.obj) <- 'DNAm'
marks <- FindAllMarkers(seurat.obj)
sel_markers <- as.vector(unlist(sapply(unique(marks$cluster), function(x)row.names(marks[marks$cluster%in%x, ])[1:2])))
sel_markers <- gsub("[[:punct:]][[:alnum:]]", "", sel_markers)
to_plot <- as.data.frame(joint.obj[['umap']]@cell.embeddings)
to_plot <- data.frame(to_plot, t(GetAssayData(joint.obj)))
to_plot <- to_plot[, c('UMAP_1',
                       'UMAP_2',
                       sel_markers)]
to_plot <- reshape2::melt(to_plot, id=c('UMAP_1', 'UMAP_2'))
to_plot$value <- as.numeric(to_plot$value)
plot_list <- list()
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
for(ab in unique(to_plot$variable)){
  sel_dat <- subset(to_plot, subset=variable==ab)
  plot <- ggplot(sel_dat, aes(x=UMAP_1, y=UMAP_2, color=value))+geom_point(size=.1, stroke=.1)+
    plot_theme_facet+scale_color_viridis(option='mako', direction=-1)+ggtitle(ab)
  plot_list[[ab]] <- plot
}
#plot_list[[8]] <- plot_list[[4]] 
#plot_list[[4]] <- ggplot()+theme_void()
png(file.path(out.folder,"UMAP_protein_methylome.png"), width = 40, height = 45, units='mm', res=300)
plots <- do.call(grid.arrange, c(plot_list, ncol = 7))
dev.off()

cts_seurat <- seurat.obj[[]]$predicted.ct
cts_sergio <- Idents(sergio.obj)
meth_data <- as.matrix(GetAssayData(seurat.obj, assay='DNAm'))
rna_data <- as.matrix(GetAssayData(sergio.obj, assay='RNA'))
cluster_mean_meth <- aggregate(t(meth_data), by=list(cts_seurat),  mean, na.rm=TRUE)
cluster_mean_rna <- aggregate(t(rna_data), by=list(cts_sergio),  mean, na.rm=TRUE)
rownames(cluster_mean_meth) <- cluster_mean_meth$Group.1
cluster_mean_meth <- data.frame(cluster_mean_meth[,- 1])
rownames(cluster_mean_rna) <- cluster_mean_rna$Group.1
cluster_mean_rna <- data.frame(cluster_mean_rna[,- 1])
cluster_mean_rna <- cluster_mean_rna[row.names(cluster_mean_meth), ]
max_cor_gene <- c()
for(cpg in colnames(cluster_mean_meth)){
  max_cor_name <- names(which.max(apply(cluster_mean_rna[, -1], 2, function(x){
    -cor(x, cluster_mean_meth[, cpg])
  })))
  max_cor <- max(apply(cluster_mean_rna[, -1], 2, function(x){
    -cor(x, cluster_mean_meth[, cpg])
  }), na.rm=TRUE)
  max_cor_gene <- c(max_cor_gene, max_cor)
  names(max_cor_gene)[length(max_cor_gene)] <- max_cor_name
}

amplicon_info <- read.table("/users/lvelten/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt",
                            sep = '\t',
                            row.names=1,
                            header=TRUE)
row.names(amplicon_info) <- amplicon_info$Type.of.amplicon
genes <- rnb.annotation2data.frame(rnb.get.annotation('genes', 'hg19'))
ens_ids <- row.names(genes[unlist(sapply(names(max_cor_gene), function(x)grep(x, genes$symbol)[1])), ])
ens_ids <- gsub("[[:punct:]][[:alnum:]]", "", ens_ids)
info <- data.frame(CpG=colnames(cluster_mean_meth),
                   Gene=names(max_cor_gene),
                   CpG_chromosome=amplicon_info[colnames(cluster_mean_meth), 'Start.hg19'],
                   CpG_start=amplicon_info[colnames(cluster_mean_meth), 'End.hg19'],
                   CpG_end=amplicon_info[colnames(cluster_mean_meth), 'CpG_Names'],
                   Gene_chromosome=genes[ens_ids, 'Chromosome'],
                   Gene_start=genes[ens_ids, 'Start'],
                   Gene_end=genes[ens_ids, 'End'],
                   Correlation=max_cor_gene)
write.csv(info, '/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Supplement/correlated_gene_information.csv')
p1 <- FeaturePlot(seurat.obj, features='AMPL131040')
p2 <- FeaturePlot(sergio.obj, features='CXCR5', reduction='MOFAUMAP')
p3 <- DimPlot(seurat.obj, group.by='predicted.ct')+scale_color_manual(values=color_map)
p4 <- DimPlot(sergio.obj, group.by='TrianaCellType', reduction='MOFAUMAP')+scale_color_manual(values=color_map)
png('/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Supplement/correlation_ampli_gene_CXCR5.png',
    height=2500, width=4000, res=300)
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()
