library(scBFA)
library(ggplot2)
library(SingleCellExperiment)
rowinfo <- read.csv('/users/lvelten/project/Methylome/analysis/missionbio/BM/Sample11_70_percent_good_performance/tsv/rowinfo.csv',
                    row.names=1)
counts <- read.table('/users/lvelten/project/Methylome/analysis/missionbio/BM/Sample11_70_percent_good_performance/tsv/Sample11_70_percent_good_performance.barcode.cell.distribution.tsv',
                     header=TRUE,
                     row.names=1)
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_line(color='black', size=.25),
                    strip.background = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')
counts <- counts[row.names(rowinfo), ]
sce <- SingleCellExperiment(assay=list(counts=t(counts)))
bfa_model <- scBFA(scData = sce, numFactors = 2)
to_plot <- as.data.frame(bfa_model$ZZ)
to_plot <- data.frame(to_plot, rowinfo)
plot <- ggplot(to_plot, aes(x=V1, y=V2, color=CellType))+geom_point()+
  plot_theme+
  scale_color_manual(values==c('naive B-cells'='#fcbd7e',
                                                                                                                            'memory B-cells'='#fc6571',
                                                                                                                            'S2-S4 cells'='#fcef7e',
                                                                                                                            'S1 cells'='#bdfc7e'))
ggsave('/users/lvelten/project/Methylome/analysis/scTAMseq_manuscript/Supplement/Sample11_BFA.pdf', plot)