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
sample <- 'GSM5935921_BM_HhaI'
protein <- 'GSM5935922_BM_HhaI_protein'
out.folder <- '~'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.title=element_blank(),
                    axis.ticks=element_line(color='black', size=.1),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.key.size = unit(2, 'mm'), #change legend key size
                    legend.key.height = unit(2, 'mm'), #change legend key height
                    legend.key.width = unit(2, 'mm'), #change legend key width
                    legend.title = element_text(color='black',size=6), #change legend title font size
                    legend.text = element_text(color='black',size=5),
                    legend.position='none',
                    plot.title=element_blank())
dat <- read.table(paste0('../../data/', sample, '.tsv.gz'),
                  header = T)
rowinfo <- read.csv(paste0("../../misc/", sample, "/tsv/rowinfo.csv"),
                    row.names = 1)
colinfo <- read.csv(paste0("../../misc/", sample, "/tsv/colinfo.csv"),
                    row.names = 1)
out.folder <- '~'
dat <- dat[row.names(rowinfo), row.names(colinfo)]
dat <- ifelse(dat>0, 1, 0)
seurat.obj <- CreateSeuratObject(t(dat),
                                 assay = "DNAm",
                                 meta.data = rowinfo)
dat <- read.table(paste0('../../data/', sample, '.tsv.gz'),
                     header = T)
ampli_info <- read.table('../../misc/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
rowinfo <- read.csv(paste0('../../misc/', sample, '/tsv/rowinfo.csv'),
                    row.names = 1)
dat <- dat[row.names(rowinfo), row.names(ampli_info)[ampli_info$Type.of.amplicon%in%'CpG.B.cell.diff']]
out.folder <- '~'
prot_data_all <- read.table(paste0('../../data/', protein, '.tsv.gz'),
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

# Obtained from https://figshare.com/ndownloader/files/28408917
sergio.obj <- readRDS('../../data/WTA_projected.rds')
Idents(sergio.obj) <- 'ct'
cells <- Idents(sergio.obj)
DefaultAssay(sergio.obj) <- 'RNA'
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
png(file.path(out.folder,"UMAP_protein_methylome.png"), width = 40, height = 45, units='mm', res=300)
plots <- do.call(grid.arrange, c(plot_list, ncol = 7))
dev.off()

sergio.obj <- readRDS('../../data/WTA_projected.rds')
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_line(color='black', size=.25),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')
ct_map <- c('HSCs & MPPs'='S1.cells',
                               'Pre-pro-B cells'='S2.cells',
                               'Pro-B cells'='S2.cells',
                               'Pre-B cells'='S2.cells',
                              'Immature B cells'='S3.S4.cells',
                                'Nonswitched memory B cells'='memory.B.cells',
                              'Class switched memory B cells'='memory.B.cells',
                               'Mature naive B cells'='naive.B.cells')
color_map <- c('NaiveB'='#fcbd7e',
               'MemoryB'='#fc3262',
               'S3S4cells'='#fcef7e',
               'S2cells'='#d7ef7e',
               'S1cells'='#bdfc7e')
Idents(sergio.obj) <- 'ct'
cts_sergio <- ct_map[as.character(Idents(sergio.obj))]
sergio.obj <- sergio.obj[, !is.na(cts_sergio)]
cts_sergio <- cts_sergio[!is.na(cts_sergio)]
meth_data <- read.csv(paste0('../../dropout_modeling/', sample, '/corrected_values.csv'),
                      row.names=1)
rna_data <- as.matrix(GetAssayData(sergio.obj, assay='RNA'))
cluster_mean_meth <- t(meth_data)
cluster_mean_rna <- aggregate(t(rna_data), by=list(cts_sergio),  mean, na.rm=TRUE)
rownames(cluster_mean_rna) <- cluster_mean_rna$Group.1
cluster_mean_rna <- data.frame(cluster_mean_rna[,- 1])
cluster_mean_rna <- cluster_mean_rna[row.names(cluster_mean_meth), ]
to_plot <- data.frame(CellType=row.names(cluster_mean_rna),
                      CXCR5=cluster_mean_rna[, 'CXCR5'],
                      AMPL131040=cluster_mean_meth[, 'AMPL131040'])
p_val <- cor.test(to_plot$CXCR5, to_plot$AMPL131040)$p.value
plot <- ggplot(to_plot, aes(x=AMPL131040, y=CXCR5, color=CellType))+geom_point()+plot_theme+
  scale_color_manual(values=color_map)+annotate(geom='text', label=paste('R²:', format(p_val, digits=3)), x=0.5, y=0.35)
ggsave(file.path(out.folder, 'CXCR5_AMPL131040_scatterplot.pdf'),
       plot,
       width=35,
       height=40,
       unit='mm')

to_plot <- data.frame(CellType=row.names(cluster_mean_rna),
                      SMARCA4=cluster_mean_rna[, 'SMARCA4'],
                      AMPL131282=cluster_mean_meth[, 'AMPL131282'])
plot <- ggplot(to_plot, aes(x=AMPL131282, y=SMARCA4, color=CellType))+geom_point()+plot_theme+
  scale_color_manual(values=color_map)
ggsave(file.path(out.folder, 'SMARCA4_AMPL131282_scatterplot.pdf'),
       plot,
       width=35,
       height=40,
       unit='mm')


cts_seurat <- ct_map[seurat.obj[[]]$predicted.ct]
Idents(sergio.obj) <- 'ct'
cts_sergio <- ct_map[as.character(Idents(sergio.obj))]
meth_data <- read.csv(paste0('../../dropout_modeling/', sample, '/corrected_values.csv'),
                      row.names=1)
rna_data <- as.matrix(GetAssayData(sergio.obj, assay='RNA'))
cluster_mean_meth <- t(meth_data)
rna_data <- as.matrix(GetAssayData(sergio.obj, assay='RNA'))
cluster_mean_rna <- aggregate(t(rna_data), by=list(cts_sergio),  mean, na.rm=TRUE)
cluster_mean_meth <- data.frame(cluster_mean_meth[,- 1])
rownames(cluster_mean_rna) <- cluster_mean_rna$Group.1
cluster_mean_rna <- data.frame(cluster_mean_rna[,- 1])
cluster_mean_rna <- cluster_mean_rna[row.names(cluster_mean_meth), ]
amplicon_info <- read.table("../../misc/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt",
                            sep = '\t',
                            row.names=1,
                            header=TRUE)
row.names(amplicon_info) <- amplicon_info$Type.of.amplicon
genes <- rnb.annotation2data.frame(rnb.get.annotation('genes', 'hg19'))
ens_ids <- row.names(genes[unlist(sapply(colnames(cluster_mean_rna)[-1], function(x)grep(x, genes$symbol)[1])), ])
ens_ids <- gsub("[[:punct:]][[:alnum:]]", "", ens_ids)
genes <- genes[ens_ids, ]
max_cor_gene <- c()
for(i in 1:length(colnames(cluster_mean_meth))){
  cpg <- colnames(cluster_mean_meth)[i]
  chr <- ampli_info[cpg, 'chr']
  str <- ampli_info[cpg, 'amplicon_start']
  end<- ampli_info[cpg, 'amplicon_end']
  stra <- genes[ens_ids, 'Strand']
  sel_genes <- genes[, 'Chromosome']==chr &
    str > genes[, 'Start']-25000&
    end < genes[, 'End']+25000
  sel_genes[is.na(sel_genes)] <- FALSE
  sel_genes_symbols <- intersect(genes[sel_genes, 'symbol'], colnames(cluster_mean_rna))
  if(length(sel_genes_symbols)==0){
    max_cor_name <- NA
    max_cor <- NA
  }else{  
    max_cor_name <- names(which.max(apply(cluster_mean_rna[, sel_genes_symbols, drop=FALSE], 2, function(x){
      -cor(x, cluster_mean_meth[, cpg])
    })))
    max_cor <- max(apply(cluster_mean_rna[, sel_genes_symbols, drop=FALSE], 2, function(x){
      -cor(x, cluster_mean_meth[, cpg])
    }), na.rm=TRUE)
  }
  max_cor_gene <- c(max_cor_gene, max_cor)
  names(max_cor_gene)[length(max_cor_gene)] <- max_cor_name
}
max_cor_gene[max_cor_gene<0] <- NA
amplicon_info <- read.table("../../misc/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt",
                            sep = '\t',
                            row.names=1,
                            header=TRUE)
row.names(amplicon_info) <- amplicon_info$Type.of.amplicon
genes <- rnb.annotation2data.frame(rnb.get.annotation('genes', 'hg19'))
ens_ids <- row.names(genes[unlist(sapply(names(max_cor_gene), function(x)grep(x, genes$symbol)[1])), ])
ens_ids <- gsub("[[:punct:]][[:alnum:]]", "", ens_ids)
info <- data.frame(CpG=colnames(cluster_mean_meth),
                   Gene=names(max_cor_gene),
                   CpG_chromosome_hg19=amplicon_info[colnames(cluster_mean_meth), 'Start.hg19'],
                   CpG_start_hg19=amplicon_info[colnames(cluster_mean_meth), 'End.hg19'],
                   CpG_end_hg19=amplicon_info[colnames(cluster_mean_meth), 'CpG_Names'],
                   Gene_chromosome_hg19=genes[ens_ids, 'Chromosome'],
                   Gene_start_hg19=genes[ens_ids, 'Start'],
                   Gene_end_hg19=genes[ens_ids, 'End'],
                   Correlation=-max_cor_gene)
write.csv(info, file.path(out.folder, 'correlated_gene_information.csv'))

sel_cts <- c('HSCs & MPPs', 'Pre-B cells', 'Pro-B cells', 'Pre-pro-B cells', 'Immature B cells',
             'Mature naive B cells', 'Class switched memory B cells', 'Nonswitched memory B cells')
sel_cells <- names(Idents(sergio.obj)[Idents(sergio.obj)%in%sel_cts])

to_plot <- as.data.frame(seurat.obj[['umap']]@cell.embeddings)
to_plot <- data.frame(to_plot, AMPL131040=t(GetAssayData(seurat.obj))[, 'AMPL131040'])
p1 <- ggplot(to_plot, aes(x=UMAP_1, y=UMAP_2, color=AMPL131040))+geom_point(size=.15, stroke=.15)+
  plot_theme_facet+scale_color_viridis(option='mako', direction=-1)+ggtitle('AMPL131040')
to_plot <- as.data.frame(sergio.obj[['Projected']]@cell.embeddings[sel_cells, ])
to_plot <- data.frame(to_plot, CXCR5=t(GetAssayData(sergio.obj))[sel_cells, 'CXCR5'])
p2 <- ggplot(to_plot, aes(x=Projected_1, y=Projected_2, color=CXCR5))+geom_point(size=.15, stroke=.15)+
  plot_theme_facet+scale_color_viridis(option='mako', direction=-1)+ggtitle('CXCR5')
png('correlation_ampli_gene_CXCR5.png',
    height=453, width=792, res=300)
grid.arrange(p1, p2, ncol=2)
dev.off()

to_plot <- as.data.frame(seurat.obj[['umap']]@cell.embeddings)
to_plot <- data.frame(to_plot, AMPL131282=t(GetAssayData(seurat.obj))[, 'AMPL131282'])
p1 <- ggplot(to_plot, aes(x=UMAP_1, y=UMAP_2, color=AMPL131282))+geom_point(size=.15, stroke=.15)+
  plot_theme_facet+scale_color_viridis(option='mako', direction=-1)+ggtitle('AMPL131282')
to_plot <- as.data.frame(sergio.obj[['Projected']]@cell.embeddings[sel_cells, ])
to_plot <- data.frame(to_plot, SMARCA4=t(GetAssayData(sergio.obj))[sel_cells, 'SMARCA4'])
p2 <- ggplot(to_plot, aes(x=Projected_1, y=Projected_2, color=SMARCA4))+geom_point(size=.15, stroke=.15)+
  plot_theme_facet+scale_color_viridis(option='mako', direction=-1)+ggtitle('SMARCA4')
png('correlation_ampli_gene_SMARCA4.png',
    height=453, width=792, res=300)
grid.arrange(p1, p2, ncol=2)
dev.off()

to_plot <- as.data.frame(sergio.obj[['Projected']]@cell.embeddings)
to_plot <- data.frame(to_plot, CellType=Idents(sergio.obj))
p2 <- ggplot(to_plot, aes(x=Projected_1, y=Projected_2, color=CellType))+geom_point(size=.15, stroke=.15)+
  plot_theme_facet+scale_color_manual(values=color_map)
png('Triana_UMAP.png',
    height=453, width=396, res=300)
p2
dev.off()

genes <- rnb.annotation2data.frame(rnb.get.annotation('genes', 'hg19'))
ens_ids <- row.names(genes[unlist(sapply(colnames(cluster_mean_rna)[-1], function(x)grep(x, genes$symbol)[1])), ])
ens_ids <- gsub("[[:punct:]][[:alnum:]]", "", ens_ids)
genes <- genes[ens_ids, ]
max_cor_gene <- c()
for(i in 1:length(colnames(cluster_mean_meth))){
  cpg <- colnames(cluster_mean_meth)[i]
  chr <- ampli_info[cpg, 'chr']
  str <- ampli_info[cpg, 'amplicon_start']
  end<- ampli_info[cpg, 'amplicon_end']
  stra <- genes[ens_ids, 'Strand']
  sel_genes <- genes[, 'Chromosome']==chr &
    str > genes[, 'Start']-25000&
    end < genes[, 'End']+25000
  sel_genes[is.na(sel_genes)] <- FALSE
  sel_genes_symbols <- intersect(genes[sel_genes, 'symbol'], colnames(cluster_mean_rna))
  if(length(sel_genes_symbols)==0){
    max_cor_name <- NA
    max_cor <- NA
  }else{  
    max_cor_name <- names(which.max(apply(cluster_mean_rna[, sel_genes_symbols, drop=FALSE], 2, function(x){
      cor(x, cluster_mean_meth[, cpg])
    })))
    max_cor <- max(apply(cluster_mean_rna[, sel_genes_symbols, drop=FALSE], 2, function(x){
      cor(x, cluster_mean_meth[, cpg])
    }), na.rm=TRUE)
  }
  max_cor_gene <- c(max_cor_gene, max_cor)
  names(max_cor_gene)[length(max_cor_gene)] <- max_cor_name
}
max_cor_gene[max_cor_gene<0] <- NA
amplicon_info <- read.table("../../misc/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt",
                            sep = '\t',
                            row.names=1,
                            header=TRUE)
row.names(amplicon_info) <- amplicon_info$Type.of.amplicon
genes <- rnb.annotation2data.frame(rnb.get.annotation('genes', 'hg19'))
ens_ids <- row.names(genes[unlist(sapply(names(max_cor_gene), function(x)grep(x, genes$symbol)[1])), ])
ens_ids <- gsub("[[:punct:]][[:alnum:]]", "", ens_ids)
info <- data.frame(CpG=colnames(cluster_mean_meth),
                   Gene=names(max_cor_gene),
                   CpG_chromosome_hg19=amplicon_info[colnames(cluster_mean_meth), 'Start.hg19'],
                   CpG_start_hg19=amplicon_info[colnames(cluster_mean_meth), 'End.hg19'],
                   CpG_end_hg19=amplicon_info[colnames(cluster_mean_meth), 'CpG_Names'],
                   Gene_chromosome_hg19=genes[ens_ids, 'Chromosome'],
                   Gene_start_hg19=genes[ens_ids, 'Start'],
                   Gene_end_hg19=genes[ens_ids, 'End'],
                   Correlation=max_cor_gene)
write.csv(info, file.path(out.folder, 'correlated_gene_information_positive.csv'))

DefaultAssay(sergio.obj) <- 'RNA'
markers <- FindMarkers(sergio.obj, ident.1=c('HSCs & MPPs', 'Pre-B cells', 'Pro-B cells', 'Pre-pro-B cells', 'Immature B cells'),
                       ident.2=c('Mature naive B cells', 'Class switched memory B cells', 'Nonswitched memory B cells'))
int_genes <- intersect(info$Gene, row.names(markers))
info[info$Gene%in%int_genes&info$Correlation>0, ]

FeaturePlot(joint.obj, features=c('AMPL130599', 'GNB1'))
p1 <- FeaturePlot(seurat.obj, features='AMPL130599')
p2 <- FeaturePlot(sergio.obj, features='GNB1', reduction='Projected', cells=colnames(sergio.obj)[Idents(sergio.obj)%in%sel_cts])
p3 <- DimPlot(seurat.obj, group.by='predicted.ct')+scale_color_manual(values=color_map)
p4 <- DimPlot(sergio.obj, group.by='ct', reduction='Projected')+scale_color_manual(values=color_map)
png('correlation_ampli_gene_GNB1_positive.png',
    height=2500, width=4000, res=300)
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

p1 <- FeaturePlot(seurat.obj, features='AMPL131282')
p2 <- FeaturePlot(sergio.obj, features='SMARCA4', reduction='Projected')
p3 <- DimPlot(seurat.obj, group.by='predicted.ct')+scale_color_manual(values=color_map)
p4 <- DimPlot(sergio.obj, group.by='ct', reduction='Projected')+scale_color_manual(values=color_map)
png('correlation_ampli_gene_SMARCA4_positive.png',
    height=2500, width=4000, res=300)
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

seurat.mono <- as.cell_data_set(joint.obj)
seurat.mono <- cluster_cells(cds = seurat.mono, reduction_method = "UMAP")
seurat.mono <- learn_graph(seurat.mono, use_partition = FALSE, close_loop=FALSE)
hscs <- names(Idents(joint.obj)[Idents(joint.obj)%in%'HSCs & MPPs'])

seurat.mono <- order_cells(seurat.mono, reduction_method = "UMAP", root_cells=hscs)

# plot trajectories colored by pseudotime
to_plot <- as.data.frame(joint.obj[['umap']]@cell.embeddings)
to_plot <- data.frame(to_plot, Pseudotime=seurat.mono@principal_graph_aux$UMAP$pseudotime)
plot <- ggplot(to_plot, aes(x=UMAP_1, y=UMAP_2, color=Pseudotime))+geom_point(size=.2, stroke=.2)+
  plot_theme+scale_color_viridis(option='inferno')
ggsave(file.path(out.folder,"UMAP_trajectory_joint.png"), width = 40, height = 45, unit='mm')
