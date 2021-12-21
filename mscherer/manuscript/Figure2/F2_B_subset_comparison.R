################### F1_F_subset_comparison.R ################### 
#' Here, we compare the subsets of memory Bcells we identify in the previous analysis
#' and compare them to bulk values of class-switch- and non-switch memory B-cells
#' 
library(ggplot2)
library(pheatmap)
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/'
cut <- 'Sample7_70_percent_good_performance'
uncut <- 'Sample6_70_percent_good_performance'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_line(color='black', size=.25),
                    strip.background = element_blank(),
#                    strip.text.x = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5),
#                    axis.text.y=element_blank(),
                    legend.position='none')
color_map <- c('naive'='#fcbd7e',
               'memory1'='#fc6571',
               'memory2'='#fc3262')

filtered.counts <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/", cut, "/tsv/", cut, ".barcode.cell.distribution.tsv"), 
                              sep="\t", 
                              header=T)
cell_metadata <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', cut, '/tsv/rowinfo.csv'),
                          row.names = 1)
amplicon.info <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/", uncut, "/tsv/selected_amplicons.tsv"),
                            row.names = 1)

selected <- ifelse(filtered.counts[row.names(cell_metadata), ]>0, 1, 0)

cell_metadata$switch <- c('Cluster1'='Cluster1',
                          'Cluster2a'='Cluster2a',
                          'Cluster2b'='Cluster2b',
                          'Cluster2c'='Cluster2c')[cell_metadata$Cluster]
# p.vals <- sapply(unique(cell_metadata$switch), function(ct1){
#   sapply(unique(cell_metadata$switch), function(ct2){
#     if(as.character(ct1)>=as.character(ct2)) return(NA)
#     unlist(apply(filtered.counts, 2, function(x, c1, c2){
#       sel_ct1 <- row.names(cell_metadata)[cell_metadata$switch%in%c1]
#       sel_ct2 <- row.names(cell_metadata)[cell_metadata$switch%in%c2]
#       sort(wilcox.test(x[sel_ct1], x[sel_ct2])$p.value)
#     }, c1=ct1, c2=ct2))
#   })
# })
# row.names(p.vals) <- colnames(p.vals) <- unique(cell_metadata$switch)
load('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/meth.data.numeric.Rdata')

#bcell_clusters <- unlist(p.vals['Cluster2b','Cluster2a&c'])
#bcell_clusters <- sort(bcell_clusters)[1:20]
bcell_clusters <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/differential/differential_CpGs_Cluster2avsCluster2b.csv')
bcell_clusters <- as.character(bcell_clusters$Amplicon)
names(bcell_clusters) <- bcell_clusters
bcell_clusters <- sort(bcell_clusters, decreasing = TRUE)[1:20]
cluster_meth <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/dropou_modeling/re_clustering/clusters/all_amplicons_clusters.csv', row.names=1)
cell_assignment <- read.table('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/MBC_assignment.txt')
meth.data.numeric <- meth.data.numeric[,as.character(cell_assignment$V5)]
more_info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
row.names(more_info) <- more_info$amplicon
joint_names <- intersect(row.names(more_info), names(bcell_clusters))
bcell_clusters <- bcell_clusters[joint_names]
cluster_meth <- cluster_meth[joint_names, ]
names(bcell_clusters) <- more_info[names(bcell_clusters), 'background.cpgs']
n.ames <- more_info[row.names(cluster_meth), 'background.cpgs']
row.names(cluster_meth) <- n.ames
cluster_meth <- cluster_meth[!is.na(row.names(cluster_meth)),]
mean_classes <- aggregate(t(meth.data.numeric[row.names(cluster_meth), ]), by=list(cell_assignment$V2), mean)
cluster_meth <- data.frame(Group.1=c('Cluster2a', 'Cluster2b'),
                           rbind(cluster_meth$Cluster2a,
                                 cluster_meth$Cluster2b))
colnames(cluster_meth) <- c('Group.1', n.ames)
to_plot <- data.frame(Type=c('Bulk', 'Bulk', 'SingleCell', 'SingleCell'),
                      rbind(mean_classes, 
                            cluster_meth))
to_plot <- reshape2::melt(to_plot, id=c('Group.1', 'Type'))
colnames(to_plot)[3:4] <- c('CpGID', 'Methylation')
to_plot$CpGID <- factor(to_plot$CpGID, levels=unique(names(sort(bcell_clusters, decreasing=TRUE))))
rename_cluster <- c('Cluster2a'='ncsMBC',
                    'Cluster2b'='csMBC',
                    'csMBC'='csMBC',
                    'ncsMBC'='ncsMBC')
to_plot$Group.1 <- rename_cluster[to_plot$Group.1]

library(viridis)
plot <- ggplot(to_plot,aes(x=Type, y=CpGID, fill=Methylation))+geom_tile()+geom_text(aes(label=format(Methylation, digits = 2)), color='white', size=2)+
  scale_fill_viridis(option='inferno', begin=1, end=0)+facet_wrap(Group.1~., strip.position = "top")+xlab('')+
  plot_theme
ggsave(file.path(plot_path, 'F2_B_subset_comparison.pdf'), plot, width=82, height=70, units='mm')
