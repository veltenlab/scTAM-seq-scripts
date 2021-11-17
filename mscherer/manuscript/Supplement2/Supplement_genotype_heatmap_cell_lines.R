.libPaths(c(.libPaths(), '/users/mscherer/conda/envs/rnbeads/lib/R/library/'))
library(snpStats)
library(VariantAnnotation)
library(ggfortify)
library(pheatmap)
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=10),
                    axis.text=element_text(color='black',size=10),
                    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    legend.key=element_rect(color='black', fill=NA))

sample <- 'Sample5_80_percent'
param <- ScanVcfParam(which = GRanges(c('chr3:194057636',
                                        'chr10:82299620',
                                       'chr1:110738296',
                                       'chr22:43910839',
                                        'chr9:137756496',
                                       'chr11:119879255',
                                        'chr4:9677287',
                                        'chr11:119879327',
                                        'chr19:1035450',
                                        'chr19:31777639',
                                        'chr16:68318995')))
vcf <- readVcf(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/vcf/',sample,'.combined.vcf.gz'),
               genome='hg19',
               param = param)
snp.mat <- genotypeToSnpMatrix(vcf)
snp.num <- as(snp.mat$genotypes,'numeric')
num.nas <- apply(snp.num,2,function(x)sum(is.na(x)))
snp.zero.nas <- snp.num[,num.nas<1]

snp.few.nas <- snp.num[,num.nas<0.05*nrow(snp.num)]
var.snps <- apply(snp.few.nas,2,var,na.rm=T)
snp.few.nas <- snp.few.nas[,var.snps>0.005]
snp.few.nas[is.na(snp.few.nas)] <- 0.0
prot.data <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'-protein-counts.tsv'),header = TRUE)
prot.data <- prot.data[(prot.data$cell_barcode%in%row.names(snp.few.nas)),]
barcodes <- unique(prot.data$cell_barcode)
cd3 <- prot.data[(prot.data$cell_barcode%in%barcodes)&(prot.data$ab_description%in%'CD3'), c('cell_barcode', 'raw')]
row.names(cd3) <- cd3$cell_barcode
cd33 <- prot.data[(prot.data$cell_barcode%in%barcodes)&(prot.data$ab_description%in%'CD33'), c('cell_barcode', 'raw')]
row.names(cd33) <- cd33$cell_barcode
prot.data <- data.frame(CD3=cd3[barcodes, 'raw'], CD33=cd33[barcodes, 'raw'])
row.names(prot.data) <- barcodes
clust <- hclust(dist(snp.few.nas),method = 'ward.D')
clusts <- cutree(clust,3)
clust <- as.factor(c('1'='K562',
           '2'='Jurkat',
           '3'='Mixed')[as.character(clusts)])
names(clust) <- names(clusts)
prot.data$cluster <- clust[row.names(prot.data)]
prot.data$CD3 <- log10(prot.data$CD3+1)
prot.data$CD33 <- log10(prot.data$CD33+1)
prot.data <- reshape2::melt(prot.data,id='cluster')
colnames(prot.data) <- c('cluster','ab_description','raw')
plot <- ggplot(prot.data,aes(x=ab_description,y=raw,fill=cluster))+geom_boxplot()+plot_theme+
  scale_fill_manual(values=c('Jurkat'='#ff9289', 'K562'='#82b7ff','Mixed'='gray'))
ggsave(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/vcf/ab_boxplots.pdf'),plot,
       width=130, height=80, unit='mm')
clust.info <- data.frame(cluster=factor(clust))
png(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/vcf/genotype_heatmap.png'), width=149, height=149, unit='mm', res=100)
pheatmap(snp.few.nas,
         show_rownames = FALSE,
         clustering_method='ward.D',
         annotation_row = clust.info,
         cutree_rows = 3,
         show_colnames = FALSE,
         legend = FALSE,
         annotation_legend = FALSE,
         annotation_colors = list(cluster=c('Jurkat'='#ff9289', 'K562'='#82b7ff','Mixed'='gray')))
dev.off()
