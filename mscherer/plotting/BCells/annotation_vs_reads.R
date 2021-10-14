library(RnBeads)
my_theme <- theme_bw()+theme(panel.grid=element_blank(),
                            text=element_text(size=18,color="black"),
                            axis.ticks=element_line(color="black"),
                            plot.title = element_text(size=18,color="black",hjust = .5),
                            axis.text = element_text(size=15,color="black"),
                            axis.text.x = element_text(angle=45,hjust=1,size=15,color="black"))


# Cell lines
ampli.info <- read.csv("/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv")
ampli.gr <- makeGRangesFromDataFrame(ampli.info[,c('chr','amplicon_start','amplicon_end')],start.field='amplicon_start',end.field='amplicon_end')
dat <- read.table("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/BCells_Sample6_70_percent/tsv/BCells_Sample6_70_percent.barcode.cell.distribution.tsv",
                     header = T)
rnb.load.annotation.from.db('ensembleRegBuildBPall')
anno <- unlist(rnb.get.annotation('ensembleRegBuildBPall'))
op <- findOverlaps(ampli.gr,anno)
to.plot <- data.frame(t(dat[,queryHits(op)]),Anno=values(anno[subjectHits(op)])$elementType)
to.plot <- reshape2::melt(to.plot,id='Anno')
plot <- ggplot(to.plot,aes(x=Anno,y=value))+geom_boxplot()+my_theme
ggsave('/users/mscherer/cluster//project/Methylome/analysis/amplicons/B_cells/R/annotation_vs_reads.pdf',plot)

cpgs <- makeGRangesFromDataFrame(rnb.annotation2data.frame(rnb.get.annotation("CpG")))
sel.anno <- anno[subjectHits(op)]
op <- findOverlaps(sel.anno,cpgs)
cpg.counts <- plyr::count(queryHits(op))
cpg.counts <- cpg.counts[order(cpg.counts$freq,decreasing=TRUE),]
to.plot <- data.frame(CpGs=cpg.counts$freq,Anno=values(sel.anno)$elementType)
plot <- ggplot(to.plot,aes(x=Anno,y=CpGs))+geom_boxplot()+my_theme

ampli.info <- read.csv("/users/lvelten/project/Methylome/infos/amplicon_info/B_cells/amplicon_info.csv")
ampli.gr <- makeGRangesFromDataFrame(ampli.info[,c('chr','amplicon_start','amplicon_end')],start.field='amplicon_start',end.field='amplicon_end')
cpgs <- makeGRangesFromDataFrame(rnb.annotation2data.frame(rnb.get.annotation("CpG")))
op <- findOverlaps(ampli.gr,anno)
sel.anno <- anno[subjectHits(op)]
op <- findOverlaps(sel.anno,cpgs)
cpg.counts <- plyr::count(queryHits(op))
cpg.counts <- cpg.counts[order(cpg.counts$freq,decreasing=TRUE),]
to.plot <- data.frame(CpGs=cpg.counts$freq,Anno=values(sel.anno)$elementType)
plot <- ggplot(to.plot,aes(x=Anno,y=CpGs))+geom_boxplot()+my_theme

