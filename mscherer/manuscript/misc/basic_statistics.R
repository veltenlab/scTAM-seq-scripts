library(ggplot2)
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_line(color='black', size=.25),
                    strip.background = element_blank(),
                    #strip.text.x = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')
sample <- 'Sample8_70_percent_good_performance'
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure1/Sample8/'

rowinfo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/', sample, '/tsv/rowinfo.csv'),
                    row.names = 1)
filtered.counts <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/", sample, "/tsv/", sample, ".barcode.cell.distribution.tsv"),
                              row.names = 1, header=T)


selected_data <- filtered.counts[row.names(rowinfo), ]
reads_per_cell <- apply(selected_data, 1, sum)
reads_per_amplicon <- apply(selected_data, 2, sum)
to_plot <- data.frame(ReadsPerCell=reads_per_cell)
plot <- ggplot(to_plot, aes(x=ReadsPerCell, y=..count..))+geom_histogram()+plot_theme+
  annotate(x=mean(reads_per_cell), y=500, color='white', geom='text', label=paste('Mean:' , format(mean(reads_per_cell, digits=0))), size=1)
ggsave(file.path(plot_path, 'reads_per_cell.pdf'), plot, width=80, height=60, unit='mm')

to_plot <- data.frame(ReadsPerAmplicon=reads_per_amplicon)
plot <- ggplot(to_plot, aes(x=ReadsPerAmplicon, y=..count..))+geom_histogram()+plot_theme+
  annotate(x=100000, y=100, color='black', geom='text', label=paste('Mean:' , format(mean(reads_per_amplicon, digits=0))), size=1)
ggsave(file.path(plot_path, 'reads_per_amplicon.pdf'), plot, width=80, height=60, unit='mm')
