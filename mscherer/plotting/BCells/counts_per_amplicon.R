############## counts_per_amplicon.R ##############
#' This file creates scatteplots for all amplicons in a given panel comparing the total counts per cell
#' versus the counts for this specifi amplicon. Additional cell clustering information has to be provided
#' in tabular form

library(dplyr)
sample <- 'BCells_Sample6_70_percent_good_performance'
dat <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/tsv/",sample,".barcode.cell.distribution.tsv"), 
                 sep="\t", 
                 header=T)
total_count_amp_cell <- apply(dat,1,sum)
log_total_count_amp_cell <- log10(total_count_amp_cell)
dat.y <- dat %>% +1 %>% log10()
pdf("/users/mscherer/cluster/project/Methylome/analysis/amplicons/B_cells/BCells_Sample6_70_percent_good_performance/amplicon_plots.pdf", width = 12, height = 5)
par(mfrow = c(2,4))
for (i in 1:ncol(dat)){
  plot(x = log_total_count_amp_cell, y= dat.y[,i], 
       ylab = "Log(Counts/ampl)",
       xlab = "log(TotalCounts)",
       pch = 20,
       cex = 0.5, cex.lab=0.7, cex.axis=0.7)
  
  title(paste(colnames(dat.y[i]),cex.main = 0.9, cex.sub = 0.8)) 
}
dev.off()
