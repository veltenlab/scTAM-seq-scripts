library(RnBeads)
library(ggsci)
library(ggfortify)
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

plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_line(color='black', size=.25),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')

load('../../data/meth.data.numeric.Rdata')
meth.data <- meth.data.numeric
fully.meth <- apply(meth.data, 1, function(x)all(x>0.75, na.rm = TRUE))
fully.unmeth <- apply(meth.data, 1, function(x)all(x<0.25, na.rm = TRUE))
rem.sites <- meth.data[!fully.meth&!fully.unmeth, ]
mean.s1 <- rowMeans(rem.sites[, grepl('S1', colnames(rem.sites))])
mean.s2 <- rowMeans(rem.sites[, grepl('S2', colnames(rem.sites))])
mean.s3 <- rowMeans(rem.sites[, grepl('S3', colnames(rem.sites))])
mean.s4 <- rowMeans(rem.sites[, grepl('S4', colnames(rem.sites))])
mean.nbc <- rowMeans(rem.sites[, grepl('NBC', colnames(rem.sites))])
mean.mbc <- rowMeans(rem.sites[, grepl('MBC', colnames(rem.sites))])
mean.gc <- rowMeans(rem.sites[, grepl('GC', colnames(rem.sites))])
all_means <- list(mean.s1,
                  mean.s2,
                  mean.s3,
                  mean.s4,
                  mean.nbc,
                  mean.mbc,
                  mean.gc)
all_differences <- as.data.frame(lapply(all_means, function(x){
  as.data.frame(lapply(all_means, function(y){
    abs(x-y)
  }))
}))
int.sites <- apply(all_differences, 1, function(x){
  any(x>0.5, na.rm=T)
})
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
                    axis.text.x=element_blank(),
                    axis.title.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    legend.title=element_blank(),
                    panel.spacing = unit(.1, "lines"))
plot <- ggplot(to_plot, aes(x="", y=value, fill=variable))+geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y", start=0)+plot_theme+scale_fill_tron()+ylab("")+xlab("")
ggsave('analyzable_CpGs_BCells.pdf',
       height=300,
       width=300,
       unit='mm')

int.sites <- names(int.sites[int.sites])
pca.obj <- prcomp(t(na.omit(meth.data)))
rowinfo <- read.table('../../misc/Code_array_samples_updated_220318.txt', header=T)
row.names(rowinfo) <- rowinfo$Sample
cols <- c('Naive.B.cell'='#fcbd7e',
          'csMemory.B.cell'='#fc6571',
          'ncsMemory.B.cell'='#fc6571',
          'PreBI'='#d7ef7e',
          'PreBII'='#d7ef7e',
          'immature.B'='#fcef7e',
          'Hematopoietic.Progenitor'='#bdfc7e')
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_line(color='black', size=.25),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')
pdf('PCA_all_CpGs.pdf',
    height=2.6,
    width=3.2)
autoplot(pca.obj, data=rowinfo, colour='Cell.type')+scale_color_manual(values=cols)+plot_theme
dev.off()

pca.red <- prcomp(t(na.omit(meth.data[int.sites, ])))
pdf('PCA_dynamic_CpGs.pdf',
    height=2.6,
    width=3.2)
autoplot(pca.red, data=rowinfo, colour='Cell.type')+scale_color_manual(values=cols)+plot_theme
dev.off()

pca.rand <- prcomp(t(na.omit(meth.data[sample(which(row.names(meth.data)%in%int.sites), 500), ])))
pdf('PCA_500_CpGs.pdf',
    height=2.6,
    width=3.2)
autoplot(pca.rand, data=rowinfo, colour='Cell.type')+scale_color_manual(values=cols)+plot_theme
dev.off()

pve.rand <- sapply(1:100, function(x){
  pca.rand <- prcomp(t(na.omit(meth.data[sample(which(row.names(meth.data)%in%int.sites), 500), ])))
  pca.rand$sdev^2 / sum(pca.rand$sdev^2)
})
pve.rand.mean <- rowMeans(pve.rand)
pve.rand.se <- apply(pve.rand, 1, sd)/sqrt(ncol(pve.rand))
pve.red <- pca.red$sdev^2 / sum(pca.red$sdev^2)
pve.all <- pca.obj$sdev^2 / sum(pca.obj$sdev^2)
to_plot <- data.frame(PC=1:length(pve.all), All=pve.all, Reduced=pve.red)
to_plot <- reshape2::melt(to_plot, id='PC')
pdf('PVE_comparison.pdf',
    height=2.6,
    width=3.36)
ggplot(to_plot, aes(x=PC, y=value, color=variable))+geom_point(size=.5)+geom_line()+plot_theme+
  ylim(0,1)+xlim(0, 65)+scale_color_manual(values=c('All'='black',
                                                    'Reduced'='#95cc5e'))
dev.off()
to_plot <- data.frame(PC=1:length(pve.all), Rand=pve.rand.mean, Rand_SE=pve.rand.se)
pdf('PVE_random.pdf',
    height=2.6,
    width=3.36)
ggplot(to_plot, aes(x=PC, y=Rand, ymax=Rand+2*Rand_SE, ymin=Rand-2*Rand_SE))+geom_errorbar(color='#709946')+geom_line(color='#709946')+plot_theme+
  ylim(0,1)+xlim(0, 65)
dev.off()

