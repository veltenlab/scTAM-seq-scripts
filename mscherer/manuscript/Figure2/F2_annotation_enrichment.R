qtlAnnotationEnrichment <- function(cpgs, background, annotation){
  ensembl.anno <- c("ctcf","distal","dnase","proximal","tfbs","tss")
  if(!(annotation%in%c("cpgislands","genes","promoters",ensembl.anno))){
    stop(paste("Invalid value for annotation, needs to be",c("cpgislands","genes","promoters",ensembl.anno)))
  }
  if(annotation%in%c("cpgislands","genes","promoters")){
    anno.type <- annotation
  }else{
    anno.type <- paste0("ensembleRegBuildBP",annotation)
    if(!anno.type%in%rnb.region.types()){
      rnb.load.annotation.from.db(anno.type)
    }
  }
  all.cpgs <- makeGRangesFromDataFrame(rnb.annotation2data.frame(rnb.get.annotation('probes450')))
  all.qtl <- all.cpgs[cpgs]
  all.input <- all.cpgs[background]
  annotation <- unlist(rnb.get.annotation(anno.type))
  tps <- length(findOverlaps(all.qtl,annotation))
  fps <- length(all.qtl)-tps
  fns <- length(findOverlaps(all.input,annotation))
  tns <- length(all.input)-fns
  gr <- fisher.test(matrix(c(tps,fns,fps,tns),2,2),alternative="greater")$p.value
  le <- fisher.test(matrix(c(tps,fns,fps,tns),2,2),alternative="less")$p.value
  or <- (tps/fps)/(fns/tns)
  return(list(enrichment=gr,depletion=le,OddsRatio=or))
}

qtlPlotAnnotationEnrichment <- function(cpgs, background){
  plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                      panel.grid=element_blank(),
                      text=element_text(color='black',size=10),
                      axis.text=element_text(color='black',size=8),
                      axis.ticks=element_line(color='black'),
                      strip.background = element_blank(),
                      strip.text.x = element_blank(),
                      legend.key=element_rect(color=NA, fill=NA),
                      axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5),
                      legend.position='none')
  all.annos <- c("cpgislands",
                 "promoters",
                 "genes",
                 "ctcf",
                 "distal",
                 "proximal",
                 "tfbs",
                 "dnase",
                 "tss")
  enr.res <- lapply(all.annos,function(anno){
      qtlAnnotationEnrichment(cpgs, background, anno)
  })
  to.plot <- data.frame(First=names(unlist(enr.res)),Second=unlist(enr.res))
  to.plot <- reshape2::melt(to.plot,id="First")
  to.plot <-data.frame(Annotation=all.annos,
                       enrichment=to.plot$value[grep("enrichment",to.plot$First)],
                       depletion=to.plot$value[grep("depletion",to.plot$First)],
                       OddsRatio=to.plot$value[grep("OddsRatio",to.plot$First)])
  plot <- ggplot(to.plot,aes(x="",y=Annotation,fill=log10(OddsRatio)))+
    geom_tile(color="black",size=ifelse(to.plot$enrichment<0.06|to.plot$depletion<0.06,5,1))+
    plot_theme+
    scale_fill_gradient2(low="dodgerblue3",mid="white",high="firebrick3")
  return(list(data=to.plot, plot=plot))
}

.libPaths(c(.libPaths(), '/users/mscherer/conda/envs/rnbeads/lib/R/library/'))
library(RnBeads)
library(LOLA)
cluster1.cpgs <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/amplicon_correlation/variable_amplicons_Cluster1.csv')
all.cpgs <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
all.cpgs <- subset(all.cpgs, subset = Type.of.amplicon=='CpG.B.cell.diff', select='background.cpgs')
res <- qtlPlotAnnotationEnrichment(cpgs=cluster1.cpgs$CpGID, all.cpgs[, 'background.cpgs'])
#all.cpgs <- row.names(rnb.annotation2data.frame(rnb.get.annotation('probes450')))
#res <- qtlPlotAnnotationEnrichment(cpgs=cluster1.cpgs$x, all.cpgs)
ggsave('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/F2_D_variable_enrichment.pdf',
       res$plot+coord_flip()+theme(axis.text.x=element_text(size=8, angle=45, vjust=1)),
       width=140.000,
       height=40,
       units='mm')

all.gr <- makeGRangesFromDataFrame(rnb.annotation2data.frame(rnb.get.annotation('probes450')))
variable.cpgs <- all.gr[cluster1.cpgs$x]
background.cpgs <- all.gr[all.cpgs[, 'background.cpgs']]
#background.cpgs <- all.gr[all.cpgs]
region.db <- loadRegionDB('/users/mscherer/cluster/project/Methylome/misc/LOLACore/hg19/')
res <- runLOLA(userSets=variable.cpgs, userUniverse=background.cpgs, regionDB = region.db)
plot <- lolaBarPlot(region.db, res, pvalCut = 1, maxTerms = 20)+plot_theme+theme(axis.text.x=element_text(size=8, angle=45, vjust=1))
ggsave('/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/F2_D_variable_enrichment_LOLA.pdf',
       plot,
       width=140.000,
       height=100,
       units='mm')

library(rtracklayer)
export.bed(variable.cpgs, '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/variable_cpgs.bed')
export.bed(background.cpgs, '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/background_cpgs.bed')

all.genes <- unlist(rnb.get.annotation('genes'))
info <- rnb.annotation2data.frame(rnb.get.annotation('genes'))
op.variable <- findOverlaps(variable.cpgs, all.genes)
variable.ids <- unlist(strsplit(info[queryHits(op.variable), 'entrezID'], ';'))
variable.ids <- variable.ids[!is.na(variable.ids)]
op.all <- findOverlaps(background.cpgs, all.genes)
all.ids <- unlist(strsplit(info[queryHits(op.all), 'entrezID'], ';'))
all.ids <- all.ids[!is.na(all.ids)]
library(GOstats)
p <- new("GOHyperGParams",
         geneIds=all.ids[1:1000],
         universeGeneIds=all.ids,
         ontology="BP",
         conditional=TRUE)
hyp <- hyperGTest(p)
ps <- probeSetSummary(hyp, 0.05, 10)

