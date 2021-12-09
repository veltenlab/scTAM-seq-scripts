##################################### Supplement_list_analyzable_CpGs.R #########################
#' With this script, we investigate which of the CpGs in the human genome can be investigated and
#' which of those are located in what kind of regulatory elements with respect to the chromatin
#' states defined in naive BCells.

.libPaths(c(.libPaths(), '/users/mscherer/conda/envs/rnbeads/lib/R/library/'))
library(RnBeads)
library(RnBeads.hg19)
library(BSgenome.Hsapiens.UCSC.hg19)
plot_path <- '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Supplement/'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=10),
                    axis.text=element_text(color='black',size=8),
                    axis.text.x=element_text(angle=90, hjust=1),
                    axis.ticks=element_line(color='black'),
                    strip.background = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')
cut.seq <- DNAString('GCGC')
genome <- BSgenome.Hsapiens.UCSC.hg19
cpgs <- rnb.get.annotation('CpG', assembly='hg19')
cpgs <- lapply(cpgs, function(x)x[seq(1, length(x), by=2)])
has.cut <- list()
total.cuts <- list()
for(chr in names(cpgs)){
    sel.cpgs <- cpgs[[chr]]
    sel.seq <- genome[[chr]]
    res <- matchPattern(cut.seq, sel.seq)
    cut.sites <- findOverlaps(sel.cpgs@ranges, res@ranges, minoverlap = 2)
    has.cut.chr <- rep(FALSE, length(sel.cpgs))
    has.cut.chr[queryHits(cut.sites)] <- TRUE
    has.cut[[chr]] <- has.cut.chr
    sel.cpgs <- resize(sel.cpgs, width = width(sel.cpgs)+300, fix='center')
    res <- matchPattern(cut.seq, sel.seq)
    cut.sites <- findOverlaps(sel.cpgs@ranges, res@ranges)
    num.cuts <- plyr::count(queryHits(cut.sites))
    num.cut.chr <- rep(0, length(sel.cpgs))
    num.cut.chr[num.cuts$x] <- num.cuts$freq
    total.cuts[[chr]] <- num.cut.chr
}
has.cut <- unlist(has.cut)
total.cuts <- unlist(total.cuts)
with.cut <- sum(has.cut)
with.single.cut <- sum(has.cut[total.cuts==1])
chrom.states <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/chromatin_states/NBCB_12_segments_intersect_hg19.bed')
chrom.states <- makeGRangesFromDataFrame(chrom.states,
                                         seqnames.field = 'V1',
                                         start.field = 'V2',
                                         end.field = 'V3',
                                         keep.extra.columns = TRUE)
info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/chromatin_states/Labels_ChromHMM.txt',
                  header=T)
map <- info$LongLabel
names(map) <- info$State
values(chrom.states)$ChromHMM <- map[values(chrom.states)$V4]
cpgs <- unlist(rnb.get.annotation('CpG', assembly='hg19'))
cpgs <- cpgs[seq(1, length(cpgs), by=2)]
all.cpgs <- length(cpgs)
op <- findOverlaps(cpgs, chrom.states)
chromatin.count.total <- plyr::count(values(chrom.states)$ChromHMM[subjectHits(op)])
op <- findOverlaps(cpgs[has.cut&total.cuts==1], chrom.states)
chromatin.count.selected <- plyr::count(values(chrom.states)$ChromHMM[subjectHits(op)])
cpgs.450 <- makeGRangesFromDataFrame(rnb.annotation2data.frame(rnb.get.annotation('probes450')))
op <- findOverlaps(cpgs.450, chrom.states)
chromatin.count.450 <- plyr::count(values(chrom.states)$ChromHMM[subjectHits(op)])
panel <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
sel.cpgs <- panel[panel$Type.of.amplicon%in%'CpG.B.cell.diff', 'background.cpgs']
sel.cpgs <- makeGRangesFromDataFrame(rnb.annotation2data.frame(rnb.get.annotation('probes450'))[sel.cpgs, ])
op <- findOverlaps(sel.cpgs, chrom.states)
chromatin.count.panel <- plyr::count(values(chrom.states)$ChromHMM[subjectHits(op)])
to_plot <- data.frame(ChromHMMState=chromatin.count.total$x,
                      AllCpGs=chromatin.count.total$freq,
                      AnalyzableCpGs=chromatin.count.selected$freq,
                      CpGs450=chromatin.count.450$freq,
                      PanelCpGs=chromatin.count.panel$freq)
to_plot <- reshape2::melt(to_plot, id='ChromHMMState')
colnames(to_plot)[2:3] <- c('Type', 'Count')
supp.labs <- c(paste("All CpGs:\n", all.cpgs), paste("Analyzable CpGs:\n", with.single.cut), paste("450k CpGs:\n", length(cpgs.450)), paste("Panel CpGs:\n", length(sel.cpgs)))
names(supp.labs) <- c("AllCpGs", "AnalyzableCpGs", 'CpGs450', 'PanelCpGs')
plot <- ggplot(to_plot, aes(x=ChromHMMState, y=Count))+geom_histogram(stat='identity')+facet_wrap(Type~., nrow=1, labeller=labeller(Type=supp.labs), scale='free')+
    plot_theme+ scale_y_continuous(labels = comma)
ggsave(file.path(plot_path, 'Supplement_analyzable_CpG.pdf'), plot)
