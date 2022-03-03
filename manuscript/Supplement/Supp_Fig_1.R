library(RnBeads)
library(ggsci)
tmp.dir <- tempdir()
tmp.file <- file.path(tmp.dir, 'RnBSet.zip')
download.file("https://rnbeads.org/materials/resource/deepTcells/import_RnBSet.zip", destfile = tmp.file)
unzip(tmp.file, exdir=tmp.dir)
rnb.set <- load.rnb.set(file.path(tmp.dir, 'import_RnBSet/'))
rem.samples <- grepl('NOMe', pheno(rnb.set)$technology)
rnb.set <- remove.samples(rnb.set, rem.samples)
meth.data <- meth(rnb.set)
covg.data <- covg(rnb.set)
min.covg <- rowMins(covg.data)
meth.data <- meth.data[min.covg>=10, ]
has.nas <- apply(meth.data, 1, function(x)any(is.na(x)))
meth.data <- meth.data[!has.nas, ]
fully.meth <- apply(meth.data, 1, function(x)all(x>0.75))
fully.unmeth <- apply(meth.data, 1, function(x)all(x<0.25))
rem.sites <- meth.data[!fully.meth&!fully.unmeth, ]
wtn.repl.diff <- rowMeans(abs(rem.sites[, c(1, 2, 3)]-rem.sites[, c(4, 5, 6)]))
btw.cell.diff <- data.frame(abs(rem.sites[, c(1, 2, 3)]-rem.sites[, c(4, 6, 1)]),
                            abs(rem.sites[, c(4, 5, 6)]-rem.sites[, c(2, 3, 1)]))
int.sites <- apply(btw.cell.diff, 1, function(x)any(x>0.5))
int.sites <- rowMaxs(as.matrix(btw.cell.diff))>wtn.repl.diff & int.sites
mean.cm <- rowMeans(rem.sites[, c(1, 4)])
mean.em <- rowMeans(rem.sites[, c(2, 5)])
mean.tn <- rowMeans(rem.sites[, c(3, 6)])
int.sites <- abs(mean.cm-mean.em) > 0.5 |
  abs(mean.cm-mean.tn) > 0.5 |
  abs(mean.tn-mean.em) > 0.5
#int.sites <- rem.sites[btw.cell.diff>wtn.repl.diff, ]
to_plot <- data.frame(FullyMethylated=sum(fully.meth),
                      FullyUnmethylated=sum(fully.unmeth),
                      Other=nrow(meth.data)-sum(fully.meth)-sum(fully.unmeth)-sum(int.sites),
                      DynamicSites=sum(int.sites))
to_plot <- reshape2::melt(to_plot)
plot_theme <- theme(panel.background = element_blank(),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_blank(),
                    strip.background = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    #axis.text.x=element_text(angle=45, hjust=1, vjust = 1),
                    axis.text.x=element_blank(),
                    axis.title.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    legend.title=element_blank(),
                    panel.spacing = unit(.1, "lines"))
plot <- ggplot(to_plot, aes(x="", y=value, fill=variable))+geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y", start=0)+plot_theme+scale_fill_tron()+ylab("")+xlab("")
ggsave('analyzable_CpGs.pdf',
       height=100,
       width=100,
       unit='mm')
