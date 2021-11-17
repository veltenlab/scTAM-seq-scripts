library(ggplot2)
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=15),
                    axis.text=element_text(color='black',size=15),
                    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    legend.key=element_rect(color='black', fill=NA))
color_map <- c('naive B-cells'='#fcbd7e',
               'memory B-cells'='#fc00bd')

filtered.counts <- read.table("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample7_70_percent_good_performance/tsv/BCells_Sample7_70_percent_good_performance.barcode.cell.distribution_with_MCL.tsv", row.names = 1, header=T)
cell_metadata <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/rowinfo_reclustering_doublet.csv',
                          row.names=1)
amplicon.info <- read.csv("/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/selected_amplicons.csv", row.names = 1)
all.amplicons <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
sel.amplicons <- row.names(all.amplicons[grepl('CpG.always.meth.B', all.amplicons$Type.of.amplicon),])
non_cut <- read.table("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample6_70_percent_good_performance/tsv/BCells_Sample6_70_percent_good_performance.barcode.cell.distribution.tsv", row.names = 1, header=T)
non_cut_doublets <- na.omit(read.csv("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample6_70_percent_good_performance/tsv/doublet_scores_DoubletDetection.csv"))
non_cut <- non_cut[non_cut_doublets[non_cut_doublets$DoubletDetectionLabel==0, 'Barcode'], ]
allelic_dropout <- apply(non_cut[, sel.amplicons], 2, function(x){
  sum(x==0)/length(x)
})
#sel.amplicons <- sel.amplicons[allelic_dropout<0.5]
allelic_dropout <- allelic_dropout[sel.amplicons]

selected <- ifelse(filtered.counts[row.names(cell_metadata), sel.amplicons]>0, 1, 0)

cell_metadata$CellType <- ifelse(cell_metadata$CellType_reclustering%in%'naive B-cells', 'naive B-cells', 'memory B-cells')
mean_methylated <- apply(selected, 2, function(x){
  sum(x>0)/length(x)
})

to_plot <- data.frame(Dropout=allelic_dropout, Methylation=mean_methylated)
plot <- ggplot(to_plot, aes(x=Dropout, y=Methylation))+geom_point()+geom_abline(slope=-1, intercept=1)

mean_methylated <- mean_methylated/(1-allelic_dropout[names(mean_methylated)])
