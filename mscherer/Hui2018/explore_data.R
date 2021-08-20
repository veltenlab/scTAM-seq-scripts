library(data.table)
all.files <- list.files('/users/mscherer/no_backup/mscherer/projects/Methylome/data/external/Hui2018/',pattern='.bed.gz')
pheno.all <- c()
pos.all <- c()
for(f in all.files){
  all.infos <- unlist(strsplit(f,'_'))
  id <- all.infos[1]
  cell.type <- all.infos[2]
  meth.data <- fread(file.path('/users/mscherer/no_backup/mscherer/projects/Methylome/data/external/Hui2018/',f))
  pos <- paste(unlist(meth.data[,1]),unlist(meth.data[,2]),unlist(meth.data[,3]),sep='_')
  pos.all <- union(pos.all,pos)
  pheno.all <- rbind(pheno.all,c(id,cell.type))
}

data.all <- matrix(NA,nrow=length(pos.all),ncol=length(all.files))
data.all <- as.data.frame(data.all)
colnames(data.all) <- pheno.all[,1]
for(f in all.files){
  all.infos <- unlist(strsplit(f,'_'))
  id <- all.infos[1]
  meth.data <- fread(file.path('/users/mscherer/no_backup/mscherer/projects/Methylome/data/external/Hui2018/',f))
  pos <- paste(unlist(meth.data[,1]),unlist(meth.data[,2]),unlist(meth.data[,3]),sep='_')
  matches.pos <- match(pos,pos.all)
  data.all[matches.pos,id] <- meth.data[,4]
}

write.csv(data.all,'/users/mscherer/cluster/project/Methylome/data/external/Hui2018/methylation_matrix.csv')

library(data.table)
library(ggplot2)
library(pheatmap)
library(ggfortify)
data.all <- fread('/users/mscherer/cluster/project/Methylome/data/external/Hui2018/methylation_matrix.csv')
pheno.all <- fread('/users/mscherer/cluster/project/Methylome/data/external/Hui2018/phenotypes.csv',header = T)
data.all <- data.all[,1:=NULL]
frac.nas.cells <- unlist(lapply(colnames(data.all),function(x){
  col <- data.all[[x]]
  sum(is.na(col))/length(col)
}))
data.all <- data.all[,(frac.nas.cells<quantile(frac.nas.cells,.9)),with=FALSE]
frac.nas <- unlist(lapply(1:nrow(data.all), function(x){
  row <- data.all[x]
  sum(is.na(row))/length(row)
}))
sel.dat <- data.all[frac.nas<0.5,]
row.names(pheno.all) <- pheno.all$SampleID
write.csv(sel.dat,'/users/mscherer/cluster/project/Methylome/data/external/Hui2018/selected_data.csv')

library(pheatmap)
sel.dat <- read.csv('/users/mscherer/cluster/project/Methylome/data/external/Hui2018/selected_data.csv')
sel.dat <- sel.dat[,-1]
pheno.all <- read.csv('/users/mscherer/cluster/project/Methylome/data/external/Hui2018/phenotypes.csv',header = T, row.names = 1)
pheatmap(as.data.frame(sel.dat),
         annotation_col = as.data.frame(pheno.all)[,c('CellType'),drop=F])
