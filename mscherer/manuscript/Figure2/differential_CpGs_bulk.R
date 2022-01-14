
load('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/meth.data.numeric.Rdata')
cell_assignment <- read.table('/users/mscherer/cluster/project/Methylome/data/external/BLUEPRINT/Renee/MBC_assignment.txt')
meth.data.numeric <- meth.data.numeric[,as.character(cell_assignment$V5)]
more_info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/CpGs.value.per.amplicon.Blood.Bone.marrow.complete.array.data.txt')
sel.cpgs <- more_info[more_info$Type.of.amplicon%in%'CpG.B.cell.diff', 'background.cpgs']
sel.amplis <- more_info[more_info$Type.of.amplicon%in%'CpG.B.cell.diff', 'amplicon']
dat_cs <- meth.data.numeric[sel.cpgs, cell_assignment$V2%in%'csMBC']
dat_ncs <- meth.data.numeric[sel.cpgs, cell_assignment$V2%in%'ncsMBC']
mean_diff <- abs(rowMeans(dat_cs)-rowMeans(dat_ncs))
p_vals <- c()
for(i in 1:nrow(dat_cs)){
  p_vals <- c(p_vals, wilcox.test(dat_cs[i, ], dat_ncs[i, ])$p.value)
}
to_write <- data.frame(Amplicon=sel.amplis,
                       CpG=sel.cpgs,
                       PValue=p_vals, 
                       MeanDiff=mean_diff)
write.csv(to_write, '/users/mscherer/cluster/project/Methylome/analysis/scTAMseq_manuscript/Figure2/differential_CpGs_bulk.csv')
