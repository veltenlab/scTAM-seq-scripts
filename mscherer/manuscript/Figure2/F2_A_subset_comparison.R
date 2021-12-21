################### F1_F_subset_comparison.R ################### 
#' Here, we compare the subsets of memory Bcells we identify in the previous analysis
#' and compare them to bulk values of class-switch- and non-switch memory B-cells
#' 
library(ggplot2)
library(pheatmap)
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/'
cut <- 'Sample7_70_percent_good_performance'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=15),
                    axis.text=element_text(color='black',size=15),
                    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=12),
                    axis.text.y=element_text(size=12),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    legend.position='none')
color_map <- c('naive'='#fcbd7e',
               'memory1'='#fc6571',
               'memory2'='#fc3262')

filtered.counts <- read.table(paste0("/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/", cut, "/tsv/", cut, ".barcode.cell.distribution.tsv"), 
                              sep="\t", 
                              header=T)
cell_metadata <- read.csv(paste0('/users/lvelten/project/Methylome/analysis/missionbio/re_sequencing/', cut, '/tsv/rowinfo.csv'),
                          row.names = 1)
amplicon.info <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/", uncut, "/tsv/selected_amplicons.tsv"),
                            row.names = 1)

selected <- ifelse(filtered.counts[row.names(cell_metadata), ]>0, 1, 0)

p.vals <- sapply(unique(cell_metadata$CellType), function(ct1){
  sapply(unique(cell_metadata$CellType), function(ct2){
    if(as.character(ct1)>=as.character(ct2)) return(NA)
    unlist(apply(filtered.counts, 2, function(x, c1, c2){
      sel_ct1 <- row.names(cell_metadata)[cell_metadata$CellType%in%c1]
      sel_ct2 <- row.names(cell_metadata)[cell_metadata$CellType%in%c2]
      sort(wilcox.test(x[sel_ct1], x[sel_ct2])$p.value)
    }, c1=ct1, c2=ct2))
  })
})
row.names(p.vals) <- colnames(p.vals) <- unique(cell_metadata$CellType)
meth_diff <- sapply(unique(cell_metadata$CellType), function(ct1){
  sapply(unique(cell_metadata$CellType), function(ct2){
    if(as.character(ct1)>=as.character(ct2)) return(NA)
    unlist(apply(filtered.counts, 2, function(x, c1, c2){
      sel_ct1 <- row.names(cell_metadata)[cell_metadata$CellType%in%c1]
      sel_ct2 <- row.names(cell_metadata)[cell_metadata$CellType%in%c2]
      meth_first <- sum(x[sel_ct1]>0)/length(sel_ct1)
      meth_second <- sum(x[sel_ct2]>0)/length(sel_ct2)
      return(c(meth_first, meth_second))
    }, c1=ct1, c2=ct2))
  })
})
row.names(meth_diff) <- colnames(meth_diff) <- unique(cell_metadata$CellType)
load('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/meth.data.numeric.Rdata')

bcell_clusters <- unlist(p.vals['memory2','memory1'])
bcell_clusters <- sort(bcell_clusters)[1:20]
bcell_clusters_diff <- t(as.data.frame(meth_diff['memory2','memory1']))
colnames(bcell_clusters_diff) <- c('memory1', 'memory2')
cell_assignment <- read.table('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/MBC_assignment.txt')
meth.data.numeric <- meth.data.numeric[,as.character(cell_assignment$V5)]
more_info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
row.names(more_info) <- more_info$amplicon
joint_names <- intersect(row.names(more_info), names(bcell_clusters))
bcell_clusters <- bcell_clusters[joint_names]
bcell_clusters_diff <- bcell_clusters_diff[joint_names,]
names(bcell_clusters) <- more_info[names(bcell_clusters), 'background.cpgs']
row.names(bcell_clusters_diff) <- more_info[row.names(bcell_clusters_diff), 'background.cpgs']
bcell_clusters_diff <- bcell_clusters_diff[!is.na(row.names(bcell_clusters_diff)),]
mean_classes <- aggregate(t(meth.data.numeric[row.names(bcell_clusters_diff), ]), by=list(cell_assignment$V2), mean)
bcell_clusters_diff <- t(bcell_clusters_diff)
bcell_clusters_diff <- data.frame(Group.1=row.names(bcell_clusters_diff), bcell_clusters_diff)
colnames(bcell_clusters_diff)[-1] <- paste0(colnames(bcell_clusters_diff)[-1])#, '/', 
                                            #more_info[match(colnames(bcell_clusters_diff)[-1], more_info$background.cpgs), 'amplicon']
colnames(mean_classes)[-1] <- paste0(colnames(mean_classes)[-1])#, '/', 
                                     #more_info[match(colnames(mean_classes)[-1], more_info$background.cpgs), 'amplicon']
to_plot <- data.frame(Type=c('Bulk', 'Bulk', 'SingleCell', 'SingleCell'),
                      rbind(mean_classes, 
                            bcell_clusters_diff))
to_plot <- reshape2::melt(to_plot, id=c('Group.1', 'Type'))
colnames(to_plot)[3:4] <- c('CpGID', 'Methylation')
to_plot$CpGID <- factor(to_plot$CpGID, levels=unique(names(sort(bcell_clusters, decreasing=TRUE))))
rename_cluster <- c('memory2'='csMBC',
                    'memory1'='ncsMBC',
                    'csMBC'='csMBC',
                    'ncsMBC'='ncsMBC')
to_plot$Group.1 <- rename_cluster[to_plot$Group.1]

library(viridis)
plot <- ggplot(to_plot,aes(x=Type, y=CpGID, fill=Methylation))+geom_tile()+geom_text(aes(label=format(Methylation, digits = 2)), color='white')+
  scale_fill_viridis(option='inferno', begin=1, end=0)+facet_wrap(Group.1~., strip.position = "top")+theme+xlab('')+
  plot_theme
ggsave(file.path(plot_path, 'F1_F_subset_comparison.pdf'), plot, width=125, height=150, units='mm')
