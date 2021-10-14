############## GC_content_vs_read_counts.R ##############
#' This file generates a plot that compares the GC content per amplicon with the read counts per amplicon across
#' all cells

library(ggplot2)
sample <- 'BCells_Sample7_70_percent_good_performance'
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
dat <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/tsv/",sample,".barcode.cell.distribution.tsv"), 
                  sep="\t", 
                  header=T)
reads.per.amplicon <- apply(dat,2,sum)
dropout.rate <- apply(dat,2,function(x){
  sum(x==0)/length(x)
})
ampli.info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
to.plot <- data.frame(DropoutRate=dropout.rate,
                      Reads=reads.per.amplicon,
                      GCContent=ampli.info[names(reads.per.amplicon),'GC_Content'])
p.val <- cor.test(to.plot$Reads,to.plot$GCContent)$p.value
plot <- ggplot(to.plot,aes(x=GCContent,y=Reads))+geom_point()+geom_smooth(method = 'lm')+
  annotate('text',label=paste('p-value: ', format(p.val,digits = 3)),x=0.32,y=750000)+theme
ggsave(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/plots/GCContent_vs_reads.pdf"),plot)

p.val <- cor.test(to.plot$DropoutRate,to.plot$GCContent)$p.value
plot <- ggplot(to.plot,aes(x=GCContent,y=DropoutRate))+geom_point()+geom_smooth(method = 'lm')+
  annotate('text',label=paste('p-value: ', format(p.val,digits = 3)),x=0.32,y=0.75)+theme
ggsave(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/plots/GCContent_vs_dropout.pdf"),plot)
