library(ggplot2)
library(data.table)
library(rtracklayer)
samples <- paste0('GSM12744',24:35)
cov.thres <- 10
orig.datas <- list.files('/users/mscherer/cluster/project/Methylome/data/external/Cabezas/',pattern = '.bed')
for(i in 1:length(samples)){
  sample <- samples[i]
  orig.data <- fread(paste0('/users/mscherer/cluster/project/Methylome/data/external/Cabezas/',orig.datas[i]))
  orig.data <- orig.data[(orig.data$V5+orig.data$V6)>cov.thres,]
  my.data <- fread(paste0('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/data/',sample,'.cov.gz'))
  my.data <- my.data[(my.data$V5+my.data$V6)>cov.thres,]
  meth.first <- as.numeric(orig.data$V4)
  meth.second <- as.numeric(my.data$V4)/100
  anno.first <- makeGRangesFromDataFrame(orig.data,seqnames.field = "V1",start.field = "V2",end.field = "V3")
  anno.second <- makeGRangesFromDataFrame(my.data,seqnames.field = "V1",start.field = "V2",end.field = "V3")
  op <- findOverlaps(anno.first,anno.second)
  to.plot <- data.frame(Original=meth.first[queryHits(op)],Replicated=meth.second[subjectHits(op)])
  cori <- round(cor(to.plot$Original,to.plot$Replicated),3)
  plot <- ggplot(to.plot,aes(x=Original,y=Replicated))+geom_point()+theme_bw()+geom_text(label=cori,x=0.1,y=0.9)+xlim(0,1)+ylim(0,1)
  ggsave(paste0('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/comparison/',sample,'_scatterplot.png'),plot)
}