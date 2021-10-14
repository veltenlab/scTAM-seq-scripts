############## ranked_cells_vs_reads.R ##############
#' This files compares the cells ranked by reads vs. the total reads per cell to determine when we start to call
#' shit as cells

library(ggplot2)
sample <- 'BCells_Sample6_70_percent_good_performance'
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
dat <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/tsv/",sample,".barcode.cell.distribution.tsv"), 
                  sep="\t", 
                  header=T)
reads.per.cell <- apply(dat,1,sum)
ranks <- nrow(dat)-rank(reads.per.cell)
working.amplicons <- apply(dat,1,function(x){
  sum(x>5)
})
to.plot <- data.frame(Rank=ranks,Reads=log10(reads.per.cell),Amplicons=working.amplicons)
#to.plot <- data.frame(Rank=ranks,Reads=reads.per.cell)
plot <- ggplot(to.plot,aes(x=Rank,y=Reads,color=Amplicons))+geom_point()+theme+ylab('log10(Number of reads per cell)')
ggsave(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/plots/read_vs_rank.pdf"),plot)
