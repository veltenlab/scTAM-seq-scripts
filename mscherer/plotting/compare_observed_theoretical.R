library(MASS)
library(ggplot2)
library(gridExtra)
sample <- 'Sample5_80_percent'
dat <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                  sep='\t',
                  header=TRUE)
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
ampli.info <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/amplicons/cell_lines/R/amplicon_info_with_RnBeads.csv')
row.names(ampli.info) <- ampli.info$AmpID_design_1898
dat <- dat[row.names(clust.file),]
dat <- dat[,row.names(ampli.info)]
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
neg.sizes.meth <- c()
neg.mu.meth <- c()
for(ampli in row.names(ampli.info)[ampli.info$Type%in%'Hha K562 high']){
  target.counts <- dat[clust.file$CellType%in%'K562',ampli]
  if(var(target.counts)<0.01) next
  total.counts <- rowSums(dat[clust.file$CellType%in%'K562',])
  to.plot <- data.frame(Target=target.counts,Total=total.counts)
  observed.plot <- ggplot(to.plot,aes(x=log10(Total),y=Target))+geom_point()+theme+ggtitle(paste('Observed:',ampli))
  params <- fitdistr(target.counts, densfun = 'negative binomial')
  neg.mu.meth <- c(neg.mu.meth,params$estimate['mu'])
  neg.sizes.meth <- c(neg.sizes.meth,params$estimate['size'])
  simulated <- rnbinom(n=length(target.counts),size=params$estimate['size'],mu=params$estimate['mu'])
  to.plot <- data.frame(Target=simulated,Total=total.counts)
  theoretical.plot <- ggplot(to.plot,aes(x=log10(Total),y=Target))+geom_point()+theme+ggtitle(paste('Theoretical:',ampli))
  to.plot <- data.frame(Theoretical=sort(simulated),Observed=sort(target.counts))
  qq.plot <- ggplot(to.plot,aes(x=Theoretical,y=Observed))+geom_point()+geom_abline(slope=1,intercept=0)+theme+ggtitle(paste('QQ-plot:',ampli))
  pdf(paste0('/users/mscherer/cluster/project/Methylome/analysis/error_mod/observed_vs_theoretical/K562_methylated/',ampli,'.pdf'))
  grid.arrange(observed.plot,theoretical.plot,qq.plot)
  dev.off()
}
neg.sizes.unmeth <- c()
neg.mu.unmeth <- c()
for(ampli in row.names(ampli.info)[ampli.info$Type%in%'Hha Jurkat high']){
  target.counts <- dat[clust.file$CellType%in%'K562',ampli]
  if(var(target.counts)<10) next
  total.counts <- rowSums(dat[clust.file$CellType%in%'K562',])
  to.plot <- data.frame(Target=target.counts,Total=total.counts)
  observed.plot <- ggplot(to.plot,aes(x=log10(Total),y=Target))+geom_point()+theme+ggtitle(paste('Observed:',ampli))
  params <- fitdistr(target.counts, densfun = 'negative binomial')
  neg.mu.unmeth <- c(neg.mu.unmeth,params$estimate['mu'])
  neg.sizes.unmeth <- c(neg.sizes.unmeth,params$estimate['size'])
  simulated <- rnbinom(n=length(target.counts),size=params$estimate['size'],mu=params$estimate['mu'])
  to.plot <- data.frame(Target=simulated,Total=total.counts)
  theoretical.plot <- ggplot(to.plot,aes(x=log10(Total),y=Target))+geom_point()+theme+ggtitle(paste('Theoretical:',ampli))
  to.plot <- data.frame(Theoretical=sort(simulated),Observed=sort(target.counts))
  qq.plot <- ggplot(to.plot,aes(x=Theoretical,y=Observed))+geom_point()+geom_abline(slope=1,intercept=0)+theme+ggtitle(paste('QQ-plot:',ampli))
  pdf(paste0('/users/mscherer/cluster/project/Methylome/analysis/error_mod/observed_vs_theoretical/K562_unmethylated/',ampli,'.pdf'))
  grid.arrange(observed.plot,theoretical.plot,qq.plot)
  dev.off()
}

neg.sizes.meth <- c()
neg.mu.meth <- c()
for(ampli in row.names(ampli.info)[ampli.info$Type%in%'Hha Jurkat high']){
  target.counts <- dat[clust.file$CellType%in%'Jurkat',ampli]
  if(var(target.counts)<0.01) next
  total.counts <- rowSums(dat[clust.file$CellType%in%'Jurkat',])
  to.plot <- data.frame(Target=target.counts,Total=total.counts)
  observed.plot <- ggplot(to.plot,aes(x=log10(Total),y=Target))+geom_point()+theme+ggtitle(paste('Observed:',ampli))
  params <- fitdistr(target.counts, densfun = 'negative binomial')
  neg.mu.meth <- c(neg.mu.meth,params$estimate['mu'])
  neg.sizes.meth <- c(neg.sizes.meth,params$estimate['size'])
  simulated <- rnbinom(n=length(target.counts),size=params$estimate['size'],mu=params$estimate['mu'])
  to.plot <- data.frame(Target=simulated,Total=total.counts)
  theoretical.plot <- ggplot(to.plot,aes(x=log10(Total),y=Target))+geom_point()+theme+ggtitle(paste('Theoretical:',ampli))
  to.plot <- data.frame(Theoretical=sort(simulated),Observed=sort(target.counts))
  qq.plot <- ggplot(to.plot,aes(x=Theoretical,y=Observed))+geom_point()+geom_abline(slope=1,intercept=0)+theme+ggtitle(paste('QQ-plot:',ampli))
  pdf(paste0('/users/mscherer/cluster/project/Methylome/analysis/error_mod/observed_vs_theoretical/Jurkat_methylated/',ampli,'.pdf'))
  grid.arrange(observed.plot,theoretical.plot,qq.plot)
  dev.off()
}
neg.sizes.unmeth <- c()
neg.mu.unmeth <- c()
for(ampli in row.names(ampli.info)[ampli.info$Type%in%'Hha K562 high']){
  target.counts <- dat[clust.file$CellType%in%'Jurkat',ampli]
  if(var(target.counts)<10) next
  total.counts <- rowSums(dat[clust.file$CellType%in%'Jurkat',])
  to.plot <- data.frame(Target=target.counts,Total=total.counts)
  observed.plot <- ggplot(to.plot,aes(x=log10(Total),y=Target))+geom_point()+theme+ggtitle(paste('Observed:',ampli))
  params <- fitdistr(target.counts, densfun = 'negative binomial')
  neg.mu.unmeth <- c(neg.mu.unmeth,params$estimate['mu'])
  neg.sizes.unmeth <- c(neg.sizes.unmeth,params$estimate['size'])
  simulated <- rnbinom(n=length(target.counts),size=params$estimate['size'],mu=params$estimate['mu'])
  to.plot <- data.frame(Target=simulated,Total=total.counts)
  theoretical.plot <- ggplot(to.plot,aes(x=log10(Total),y=Target))+geom_point()+theme+ggtitle(paste('Theoretical:',ampli))
  to.plot <- data.frame(Theoretical=sort(simulated),Observed=sort(target.counts))
  qq.plot <- ggplot(to.plot,aes(x=Theoretical,y=Observed))+geom_point()+geom_abline(slope=1,intercept=0)+theme+ggtitle(paste('QQ-plot:',ampli))
  pdf(paste0('/users/mscherer/cluster/project/Methylome/analysis/error_mod/observed_vs_theoretical/Jurkat_unmethylated/',ampli,'.pdf'))
  grid.arrange(observed.plot,theoretical.plot,qq.plot)
  dev.off()
}

test.dat <- dat[,ampli]
my.fun <- function(x, pi, mu1, mu2, theta1, theta2){
  mu1 <- exp(mu1)
  mu2 <- exp(mu2)
  theta1 <- exp(theta1)
  theta2 <- exp(theta2)
  if(pi>0.5){
    l <- dnbinom(x,size=theta1,mu=mu1)
  }else{
    l <- dnbinom(x,size=theta2,mu=mu2)
  }
  print(pi)
  print(mu1)
  print(mu2)
  print(theta1)
  print(theta2)
  print(sum(l))
  return(sum(log(l+1)))
}
fitdistr(test.dat, densfun = my.fun, start = list(pi=0,mu1=2,mu2=1,theta1=1,theta2=2), lower=c(0,1,1,1,1), upper=c(1,50,50,50,50))
