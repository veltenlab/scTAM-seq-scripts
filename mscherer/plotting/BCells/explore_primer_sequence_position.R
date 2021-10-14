#################### explore_primer_sequence_position.R ##########################
#' This files compares the sequence composition of the primers for amplicons with the dropout rate observed

library(BSgenome.Hsapiens.UCSC.hg19)
library(RnBeads)
sample <- 'BCells_Sample6_70_percent_good_performance'
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
genome <- BSgenome.Hsapiens.UCSC.hg19
dat <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/tsv/",sample,".barcode.cell.distribution.tsv"), 
                  sep="\t", 
                  header=T)
ampli.info <- read.table("/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.tsv")
ampli.gr <- makeGRangesFromDataFrame(ampli.info[,c('chr','amplicon_start','amplicon_end')],start.field='amplicon_start',end.field='amplicon_end')
dropout.per.amplicon <- apply(dat,2,function(x){
  sum(x==0)/length(x)
})

#Nucleotides
seq.content <- list()
for(i in 1:nrow(ampli.info)){
  sel.seq <- DNAString(as.character(ampli.info[i,'fwd_seq']))
  freq <- alphabetFrequency(sel.seq)
  seq.content[[row.names(ampli.info)[i]]] <- as.list(freq/length(sel.seq))
}
seq.content <- do.call(rbind.data.frame,seq.content)
to.plot <- data.frame(seq.content,dropout.per.amplicon)
plot.path <- paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/plots/sequence_content_fwd_primers/")
for(b in c('A','C','G','T')){
  cori <- cor.test(to.plot[,b],to.plot$dropout.per.amplicon)$p.value
  plot <- ggplot(to.plot,aes_string(y='dropout.per.amplicon',x=b))+geom_point()+geom_smooth(method='lm')+
    theme+annotate('text', x=0.1,y=0.4,label=format(cori,digits=3))
  ggsave(paste0(plot.path,'/',b,'.pdf'),plot)
}

#Diucleotides
seq.content <- list()
for(i in 1:nrow(ampli.info)){
  sel.seq <- DNAString(as.character(ampli.info[i,'fwd_seq']))
  freq <- dinucleotideFrequency(sel.seq)
  seq.content[[row.names(ampli.info)[i]]] <- as.list(freq/length(sel.seq))
}
seq.content <- do.call(rbind.data.frame,seq.content)
to.plot <- data.frame(seq.content,dropout.per.amplicon)
plot.path <- paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/plots/sequence_content_fwd_primers/")
for(b in colnames(seq.content)){
  cori <- cor.test(to.plot[,b],to.plot$dropout.per.amplicon)$p.value
  plot <- ggplot(to.plot,aes_string(y='dropout.per.amplicon',x=b))+geom_point()+geom_smooth(method='lm')+
    theme+annotate('text', x=0.1,y=0.4,label=format(cori,digits=3))
  ggsave(paste0(plot.path,'/',b,'.pdf'),plot)
}

seq.content <- list()
for(i in 1:nrow(ampli.info)){
  sel.seq <- DNAString(as.character(ampli.info[i,'rev_seq']))
  freq <- alphabetFrequency(sel.seq)
  seq.content[[row.names(ampli.info)[i]]] <- as.list(freq/length(sel.seq))
}
seq.content <- do.call(rbind.data.frame,seq.content)
to.plot <- data.frame(seq.content,dropout.per.amplicon)
plot.path <- paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/plots/sequence_content_bwd_primers/")
for(b in c('A','C','G','T')){
  cori <- cor.test(to.plot[,b],to.plot$dropout.per.amplicon)$p.value
  plot <- ggplot(to.plot,aes_string(y='dropout.per.amplicon',x=b))+geom_point()+geom_smooth(method='lm')+
    theme+annotate('text', x=0.1,y=0.4,label=format(cori,digits=3))
  ggsave(paste0(plot.path,'/',b,'.pdf'),plot)
}

#Diucleotides
seq.content <- list()
for(i in 1:nrow(ampli.info)){
  sel.seq <- DNAString(as.character(ampli.info[i,'rev_seq']))
  freq <- dinucleotideFrequency(sel.seq)
  seq.content[[row.names(ampli.info)[i]]] <- as.list(freq/length(sel.seq))
}
seq.content <- do.call(rbind.data.frame,seq.content)
to.plot <- data.frame(seq.content,dropout.per.amplicon)
plot.path <- paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/plots/sequence_content_bwd_primers/")
for(b in colnames(seq.content)){
  cori <- cor.test(to.plot[,b],to.plot$dropout.per.amplicon)$p.value
  plot <- ggplot(to.plot,aes_string(y='dropout.per.amplicon',x=b))+geom_point()+geom_smooth(method='lm')+
    theme+annotate('text', x=0.1,y=0.4,label=format(cori,digits=3))
  ggsave(paste0(plot.path,'/',b,'.pdf'),plot)
}

last.base <- list()
for(i in 1:nrow(ampli.info)){
  sel.seq <- as.character(ampli.info[i,'fwd_seq'])
  res <- substr(sel.seq,nchar(sel.seq),nchar(sel.seq))
  last.base[[row.names(ampli.info)[i]]] <- res
}
last.base <- do.call(rbind.data.frame,last.base)
to.plot <- data.frame(LastBase=unlist(last.base),Dropout=dropout.per.amplicon)
plot.path <- paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/plots/last_base_fwd_primer/")
plot <- ggplot(to.plot,aes(x=LastBase,y=Dropout))+geom_violin()+theme
ggsave(file.path(plot.path, 'last_base_fwd.pdf'),plot)

last.two.bases <- list()
for(i in 1:nrow(ampli.info)){
  sel.seq <- as.character(ampli.info[i,'fwd_seq'])
  res <- substr(sel.seq,nchar(sel.seq)-1,nchar(sel.seq))
  last.two.bases[[row.names(ampli.info)[i]]] <- res
}
last.two.bases <- do.call(rbind.data.frame,last.two.bases)
to.plot <- data.frame(LastBase=unlist(last.two.bases),Dropout=dropout.per.amplicon)
plot.path <- paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/plots/last_base_fwd_primer/")
plot <- ggplot(to.plot,aes(x=LastBase,y=Dropout))+geom_violin()+theme
ggsave(file.path(plot.path, 'last_two_bases_fwd.pdf'),plot)

last.three.bases <- list()
for(i in 1:nrow(ampli.info)){
  sel.seq <- as.character(ampli.info[i,'fwd_seq'])
  res <- substr(sel.seq,nchar(sel.seq)-2,nchar(sel.seq))
  last.three.bases[[row.names(ampli.info)[i]]] <- res
}
last.three.bases <- do.call(rbind.data.frame,last.three.bases)
to.plot <- data.frame(LastBase=unlist(last.three.bases),Dropout=dropout.per.amplicon)
plot.path <- paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/plots/last_base_fwd_primer/")
plot <- ggplot(to.plot,aes(x=LastBase,y=Dropout))+geom_violin()+theme
ggsave(file.path(plot.path, 'last_three_bases_fwd.pdf'),plot)

last.base <- list()
for(i in 1:nrow(ampli.info)){
  sel.seq <- as.character(ampli.info[i,'rev_seq'])
  res <- substr(sel.seq,1,1)
  last.base[[row.names(ampli.info)[i]]] <- res
}
last.base <- do.call(rbind.data.frame,last.base)
to.plot <- data.frame(LastBase=unlist(last.base),Dropout=dropout.per.amplicon)
plot.path <- paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/plots/last_base_bwd_primer/")
plot <- ggplot(to.plot,aes(x=LastBase,y=Dropout))+geom_violin()+theme
ggsave(file.path(plot.path, 'last_base_bwd.pdf'),plot)


last.two.bases <- list()
for(i in 1:nrow(ampli.info)){
  sel.seq <- as.character(ampli.info[i,'rev_seq'])
  res <- substr(sel.seq,1,2)
  last.two.bases[[row.names(ampli.info)[i]]] <- res
}
last.two.bases <- do.call(rbind.data.frame,last.two.bases)
to.plot <- data.frame(LastBase=unlist(last.two.bases),Dropout=dropout.per.amplicon)
plot.path <- paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/plots/last_base_bwd_primer/")
plot <- ggplot(to.plot,aes(x=LastBase,y=Dropout))+geom_violin()+theme
ggsave(file.path(plot.path, 'last_two_bases_bwd.pdf'),plot)

last.three.bases <- list()
for(i in 1:nrow(ampli.info)){
  sel.seq <- as.character(ampli.info[i,'rev_seq'])
  res <- substr(sel.seq,1,3)
  last.three.bases[[row.names(ampli.info)[i]]] <- res
}
last.three.bases <- do.call(rbind.data.frame,last.three.bases)
to.plot <- data.frame(LastBase=unlist(last.three.bases),Dropout=dropout.per.amplicon)
plot.path <- paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/plots/last_base_bwd_primer/")
plot <- ggplot(to.plot,aes(x=LastBase,y=Dropout))+geom_violin()+theme
ggsave(file.path(plot.path, 'last_three_bases_bwd.pdf'),plot)
