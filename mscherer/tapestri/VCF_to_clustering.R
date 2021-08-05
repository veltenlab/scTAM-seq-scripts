.libPaths(c(.libPaths(),'/users/mscherer/R'))
library(snpStats)
library(VariantAnnotation)
library(ggfortify)
library(pheatmap)
sample <- 'Sample5_80_percent'
param <- ScanVcfParam(which = GRanges(c('chr3:194057636',
#                                        'chr10:82299620',
#                                       'chr1:110738296',
#                                       'chr22:43910839',
                                        'chr9:137756496',
#                                       'chr11:119879255',
#                                        'chr4:9677287',
#                                        'chr11:119879327',
#                                        'chr12:133347014',
#                                        'chr19:1035450',
 #                                       'chr19:31777639',
                                      'chr16:68318995')))#,
   #                                     'chr9:137756645')))
vcf <- readVcf(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/vcf/',sample,'.combined.vcf.gz'),
               genome='hg19',
               param = param)
snp.mat <- genotypeToSnpMatrix(vcf)
snp.num <- as(snp.mat$genotypes,'numeric')
num.nas <- apply(snp.num,2,function(x)sum(is.na(x)))
snp.zero.nas <- snp.num[,num.nas<1000]
#pca.obj <- prcomp(t(snp.zero.nas))
#autoplot(pca.obj)
#pheatmap(snp.num)
pheatmap(snp.zero.nas)

snp.few.nas <- snp.num[,num.nas<0.05*nrow(snp.num)]
var.snps <- apply(snp.few.nas,2,var,na.rm=T)
snp.few.nas <- snp.few.nas[,var.snps>0.005]
snp.few.nas[is.na(snp.few.nas)] <- 0.0
prot.data <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'-protein-counts.tsv'),
                        header = TRUE,
                        sep='\t')
cd3 <- rep('NA',length(unique(prot.data$cell_barcode)))
cd33 <- cd3
names(cd3) <- names(cd33) <- unique(prot.data$cell_barcode)
cd3[prot.data$cell_barcode[prot.data$ab_description%in%'CD3']] <- prot.data$raw[prot.data$ab_description%in%'CD3']
cd33[prot.data$cell_barcode[prot.data$ab_description%in%'CD33']] <- prot.data$raw[prot.data$ab_description%in%'CD33']
prot.data.new <- data.frame(CD3=as.numeric(cd3),CD33=as.numeric(cd33))
row.names(prot.data.new) <- unique(prot.data$cell_barcode)
prot.data <- prot.data.new
clust <- hclust(dist(snp.few.nas),method = 'ward.D')
clust <- cutree(clust,3)
#row.names(prot.data) <- paste0(row.names(prot.data),'-1')
matchi <- match(row.names(prot.data),names(clust))
prot.data$cluster <- clust[matchi]
prot.data <- prot.data[!is.na(prot.data$cluster),]
prot.data$CD3 <- log10(prot.data$CD3+1)
prot.data$CD33 <- log10(prot.data$CD33+1)
prot.data <- reshape2::melt(prot.data,id='cluster')
colnames(prot.data) <- c('cluster','ab_description','raw')
plot <- ggplot(prot.data,aes(x=ab_description,y=raw,fill=factor(cluster)))+geom_boxplot()+theme_minimal()+
  scale_fill_manual(values=c('1'='#00ba38','2'='#f8766d','3'='#619cff'))
ggsave(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/ab_boxplots.png'),plot)
clust.info <- as.data.frame(factor(clust))
pheatmap(snp.few.nas,
         show_rownames = FALSE,
         clustering_method='ward.D',
         annotation_row = clust.info)

#clust <- 4 - clust
cluster.mean.first <- mean(prot.data$raw[prot.data$ab_description=='CD3'&prot.data$cluster==1])
cluster.mean.second <- mean(prot.data$raw[prot.data$ab_description=='CD3'&prot.data$cluster==2])
is.jurkat <- ifelse(cluster.mean.first>cluster.mean.second,1,2)
is.k562 <- ifelse(cluster.mean.first>cluster.mean.second,2,1)

is.jurkat <- 2
is.k562 <- 1
cell.assignment <- rep('Mixed',length(clust))
cell.assignment[clust==is.jurkat] <- 'Jurkat'
cell.assignment[clust==is.k562] <- 'K562'
to.write <- data.frame(Barcode=names(clust),Cluster=clust,CellType=cell.assignment)
write.table(to.write,paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'),sep='\t',quote = F)
