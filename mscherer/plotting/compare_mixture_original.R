############## compare_mixture_original.R ##############
#' Through this script, we validate whether the output of the autoencoder makes sense. To that end, we create plots
#' (scatterplots, etc.) to compare the reads counts in the original data with the output of the autoencoder

library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggplotify)
library(viridis)
sel.amplis <- FALSE
drop.cells <- FALSE
sample <- 'Sample5_80_percent'
plot.path <- paste0('/users/mscherer/cluster/project/Methylome/analysis/error_mod/compare_read_counts/',sample)
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
theme_small <- theme(panel.background = element_rect(color='black',fill='white'),
                   panel.grid=element_blank(),
                   text=element_text(color='black',size=3),
                   axis.ticks=element_line(color='black',size=1),
                   axis.ticks.length=unit(.1, "mm"),
                   legend.position='none',
                   axis.title=element_blank())

meth.data <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/methylation_autoencoder/mixture_prob_dca.csv'),
                      header = TRUE,
                      row.names = 1)
input <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                    sep='\t',
                    header=TRUE)

clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
if(drop.cells){
  rem.cells <- clust.file$Barcode[clust.file$CellType=='Mixed']
  input <- input[-which(row.names(input)%in%rem.cells),]
}
ampli.info <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/amplicons/cell_lines/R/amplicon_info_with_RnBeads.csv')
row.names(ampli.info) <- ampli.info$AmpID_design_1898
if(sel.amplis){
  rel.info <- read.csv('/users/mscherer/cluster/project/Methylome/infos/cell_lines/analysis_Agostina_Apr21/Summary_methylation_values_amplicons_K562_Jurkat.csv')
  high.conf <- rel.info$AmpID_design_1898[grepl('Jurkat.low.K562.high|Jurkat.high.K562.low',rel.info$Label)]
  ampli.info <- ampli.info[high.conf,]
  input <- input[,high.conf]
}
colnames(meth.data) <- colnames(input)
row.names(meth.data) <- row.names(input)
cmd <- paste("mkdir",plot.path)
system(cmd)
clust.file <- clust.file[row.names(meth.data),]
all.plots <- list()
for(ampli in colnames(meth.data)){
  type <- ampli.info[ampli,'Type']
  type <- gsub(' ','_',type)
  to.plot <- data.frame(Pi=meth.data[,ampli],ReadCounts=input[,ampli], CellType=clust.file$CellType)
  plot <- ggplot(to.plot,aes(x=Pi,y=log10(ReadCounts+1)))+geom_point(aes(color=CellType))+geom_smooth(method='lm')+theme+
    scale_color_manual(values=c('Jurkat'='#ff9289', 'K562'='#82b7ff', 'Mixed'='grey'))
  ampli.id <- paste0(ampli,'_',type)
  na <- paste0(plot.path,'/',ampli.id,'.pdf')
  ggsave(na, plot)
  if(type=='Hha_K562_high'){
    theme_used <- theme_small + theme(panel.border=element_rect(color='#82b7ff',fill=NA))
  }else if(type=='Hha_Jurkat_high'){
    theme_used <- theme_small + theme(panel.border=element_rect(color='#ff9289',fill=NA))
  }else{
    theme_used <- theme_small
  }
  plot <- ggplot(to.plot,aes(x=Pi,y=log10(ReadCounts+1)))+geom_point(aes(color=CellType),size=.1)+geom_smooth(method='lm',size=.5)+theme_used+
    scale_color_manual(values=c('Jurkat'='#ff9289', 'K562'='#82b7ff', 'Mixed'='grey'))+xlim(0,1)+ggtitle(ampli.id)

  all.plots[[ampli]] <- plot
}
res <- grid.arrange(grobs = all.plots[1:100],ncol=10)
ggsave(file.path(plot.path,'combined_plot_1.png'),plot(res))
res <- grid.arrange(grobs = all.plots[101:199],ncol=10)
ggsave(file.path(plot.path,'combined_plot_2.png'),plot(res))
