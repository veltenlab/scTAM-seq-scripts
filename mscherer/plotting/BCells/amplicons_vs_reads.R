############## amplicons_vs_reads.R ##############
#' This file generates a plot comparing the number of 'working' amplicons, defined as those havin more than 'cut.off'
#' reads, with the overall number of reads per cell

library(ggplot2)
sample <- 'BCells_Sample7_70_percent_good_performance'
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
dat <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/tsv/",sample,".barcode.cell.distribution.tsv"), 
                  sep="\t", 
                  header=T)
reads.per.cell <- apply(dat,1,sum)
working.amplicons <- apply(dat,1,function(x){
  sum(x>5)
})
ranks <- nrow(dat)-rank(reads.per.cell)
to.plot <- data.frame(Reads=reads.per.cell, Amplicons=working.amplicons,CellRank=ranks)
plot <- ggplot(to.plot,aes(x=Amplicons,y=Reads,color=CellRank))+geom_point()+theme
ggsave(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/plots/amplicons_vs_reads.pdf"),plot)
