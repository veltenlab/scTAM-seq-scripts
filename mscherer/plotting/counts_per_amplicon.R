############## counts_per_amplicon.R ##############
#' This file creates scatteplots for all amplicons in a given panel comparing the total counts per cell
#' versus the counts for this specifi amplicon. Additional cell clustering information has to be provided
#' in tabular form

library(dplyr)
sample <- 'Sample5_80_percent'
dat <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/tsv/",sample,".barcode.cell.distribution.tsv"), 
                 sep="\t", 
                 header=T)
clust <- read.table(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/tsv/cluster_assignment.tsv"), 
                     sep="\t", 
                     header=T)
clust <- clust[row.names(dat),]
clust$Color <- rep('steelblue3',nrow(clust))
clust$Color[clust$CellType%in%'K562'] <- 'indianred2'
clust$Color[clust$CellType%in%'Mixed'] <- 'gray'
total_count_amp_cell <- apply(dat,1,sum)
log_total_count_amp_cell <- log10(total_count_amp_cell)
dat.y <- dat %>% +1 %>% log10()
pdf("/users/mscherer/cluster/project/Methylome/analysis/amplicons/cell_lines/Sample5_80_percent/total_read_amplicons.pdf", width = 12, height = 5)
par(mfrow = c(2,4))
for (i in 1:ncol(dat)){
  plot(x = log_total_count_amp_cell, y= dat.y[,i], 
       ylab = "Log(Counts/ampl)",
       xlab = "log(TotalCounts)",
       col = as.character(clust$Color),
       pch = 20,
       cex = 0.5, cex.lab=0.7, cex.axis=0.7)
  
  title(paste(colnames(dat.y[i]),cex.main = 0.9, cex.sub = 0.8)) 
  
  lines(lowess(x = log_total_count_amp_cell[clust$Color == "steelblue3"],  y = dat.y[clust$Color == "steelblue3", i]), col= "blue")
  lines(lowess(x = log_total_count_amp_cell[clust$Color == "indianred2"],  y = dat.y[clust$Color == "indianred2", i]), col= "red")
  
  legend_vector <- c("Jurkat", "K562")
  legend("bottomright",inset=.02,box.col = NA, cex=0.75,legend_vector,fill = c("steelblue3","indianred2"))
}
dev.off()
