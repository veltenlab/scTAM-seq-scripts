library(ggplot2)
sample <- 'BCells_Sample6_70_percent_good_performance'
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
dat <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/tsv/",sample,".barcode.cell.distribution.tsv"), 
                  sep="\t", 
                  header=T)
bulk.data <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
row.names(bulk.data) <- bulk.data$amplicon
mean.counts <- apply(dat,2,mean)
plot.path <- paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/plots/bulk/')
for(ct in c('S1.mean','S2.mean','S3.mean','S4.mean','NBC.mean','GC.mean','MBC.mean','PB.mean')){
  sel.dat <- bulk.data[names(mean.counts),ct]
  type <- bulk.data[names(mean.counts),'Type.of.amplicon']
  to.plot <- data.frame(ReadCounts=mean.counts,Bulk=sel.dat,Type=type)
  cori <- cor.test(to.plot$ReadCounts,to.plot$Bulk,method='spearman')$p.value
  plot <- ggplot(to.plot,aes(x=sel.dat,y=ReadCounts,color=Type))+geom_point()+geom_smooth(method='lm')+theme+
    annotate('text',label=format(cori,digits=3),x=0.25,y=150)
  ggsave(paste0(plot.path,ct,'.pdf'),plot)
}

mean.bulk <- apply(bulk.data[,c('S1.mean','S2.mean','S3.mean','S4.mean','NBC.mean','GC.mean','MBC.mean','PB.mean')],1,mean)
dropout.rate <- apply(dat,2,function(x){
  sum(x==0)/length(x)
})
to.plot <- data.frame(Dropout=dropout.rate[names(mean.bulk)],MeanBulk=mean.bulk)
cori <- cor.test(to.plot$Dropout,to.plot$MeanBulk,method='spearman')$p.value
plot <- ggplot(to.plot,aes(x=MeanBulk,y=Dropout))+geom_point()+geom_smooth(method='lm')+theme+
  annotate('text',label=format(cori,digits=3),x=0.25,y=0.25)
ggsave(paste0(plot.path,'mean_bulk_vs_dropout.pdf'),plot)
