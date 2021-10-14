############## plot_Y_chromosome_reads.R ##############
#' This file investigates the reads on the Y chromosome to discern male from female cells

library(ggplot2)
sample <- 'BCells_Sample7_70_percent_good_performance'
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
dat <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/tsv/",sample,".barcode.cell.distribution.tsv"), 
                  sep="\t", 
                  header=T)
ampli.info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.tsv')
sel.dat <- dat[,ampli.info$chr%in%'chrY']
#AMPL202874 was not wo
sel.dat <- sel.dat[,-3]
is.male <- apply(sel.dat,1,function(x){
  any(x>1)
})
sex <- rep('female',nrow(dat))
sex[is.male] <- 'male'
to.write <- data.frame(sex)
row.names(to.write) <- row.names(dat)
write.csv(to.write,'/users/mscherer/cluster/project/Methylome/infos/BCells/cell_sexes_Sample7.csv')
