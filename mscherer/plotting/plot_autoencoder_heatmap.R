############## plot_autoencoder_heatmap.R ##############
#' With this script, we use the output of the methylation autoencoder, i.e., the mixture parameter of the forground and background read counts
#' distributions and plot it as a heatmap. The original read counts also have to be specified, as well as information about the amplicons.
#' An assignment of the cells to cell types also has to be specified and is used for cluster annotation.
#' Optionally, the matrix can be binarized to have an even cleaner picture.

library(pheatmap)
library(viridis)
sel.amplis <- FALSE
drop.cells <- FALSE
sample <- 'Sample5_80_percent'
meth.data <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/methylation_autoencoder/mixture_prob_dca_test.csv'),
                       header = TRUE,
                       row.names = 1)
input <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                   sep='\t',
                   header=TRUE)
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
if(drop.cells){
  rem.cells <- clust.file$Barcode[clust.file$CellType=='Mixed']
  input <- input[-which(row.names(input)%in%rem.cells),]
}
ampli.info <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/amplicons/cell_lines/R/amplicon_info_with_RnBeads.csv')
row.names(ampli.info) <- ampli.info$AmpID_design_1898
if(sel.amplis){
  rel.info <- read.csv('/users/mscherer/cluster/project/Methylome/infos/cell_lines/analysis_Agostina_Apr21/Summary_methylation_values_amplicons_K562_Jurkat.csv')
  high.conf <- rel.info$AmpID_design_1898[grepl('Jurkat.low.K562.high|Jurkat.high.K562.low',rel.info$Label)]
  ampli.info <- ampli.info[high.conf,]
  input <- input[,high.conf]
}
colnames(meth.data) <- colnames(input)
row.names(meth.data) <- row.names(input)
clust.file <- clust.file[row.names(meth.data),]
row.names(clust.file) <- row.names(meth.data)
ampli.info$Type[grepl('Aci',ampli.info$Type)] <- 'Aci'
ampli.info$Type <- as.factor(ampli.info$Type)
meth.data <- meth.data[,row.names(ampli.info)]
meth.data <- meth.data[,order(ampli.info$Type)]
ampli.info <- ampli.info[order(ampli.info$Type),]
for(cat in unique(ampli.info$Type)){
  sel.cols <- ampli.info$Type%in%cat
  sel.means <- apply(meth.data[,sel.cols],2,mean,na.rm=TRUE)
  meth.data[,sel.cols] <- meth.data[,sel.cols][,order(sel.means)]
}
pheatmap(meth.data,
         annotation_row = clust.file[,'CellType',drop=FALSE],
         color=rev(inferno(50)),
         annotation_col = ampli.info[,c('Type','CG'),drop=FALSE],
         show_rownames=FALSE,show_colnames=FALSE,
         clustering_method = 'ward.D',
         cluster_cols = FALSE,
         annotation_colors=list(CellType=c('Jurkat'='#ff9289', 'K562'='#82b7ff', 'Mixed'='grey'),
                                Type=c('Aci'='grey','Hha Jurkat high'='#b3645eff','Hha K562 high'='#597daeff')))

bin.map <- as.data.frame(apply(meth.data,2,round))
row.names(bin.map) <- row.names(meth.data)
colnames(bin.map) <- colnames(meth.data)
pheatmap(bin.map,
         color=rev(inferno(50)),
         annotation_row = clust.file[,'CellType',drop=FALSE],
         annotation_col = ampli.info[,c('Type','CG'),drop=FALSE],
         show_rownames=FALSE,show_colnames=FALSE,
         clustering_method = 'ward.D',
         cluster_cols = FALSE,
         annotation_colors=list(CellType=c('Jurkat'='#ff9289', 'K562'='#82b7ff', 'Mixed'='grey'),
                                Type=c('Aci'='grey','Hha Jurkat high'='#b3645eff','Hha K562 high'='#597daeff')))

library(caret)
meth.data <- apply(meth.data,2,function(x){
  x[x==0] <- 1e-5
  x[x==1] <- 0.99999
  log2(x/(1-x))
})
for.reg <- data.frame(meth.data,CellType=as.factor(ifelse(clust.file$CellType=='Mixed',1,0)))
train_control <- trainControl(method = "cv", number = 10)
model <- train(CellType ~ .,
               data = for.reg,
               trControl = train_control,
               method = "glm",
               family = binomial())
to.plot <- data.frame(for.reg,CT=clust.file$CellType)
ggplot(to.plot,aes(x=CT,y=AMPL130902))+geom_violin()+theme_bw()
sum <- summary(model)
p.vals <- sum$coefficients[,'Pr(>|z|)']
ampli.sig <- names(p.vals[p.vals<0.001])
types <- as.character(ampli.info[,'Type'])
types[match(ampli.sig,row.names(ampli.info))] <- 'significant'
ampli.info$Type <- types
pheatmap(meth.data,
         annotation_row = clust.file[,'CellType',drop=FALSE],
         color=rev(inferno(50)),
         annotation_col = ampli.info[,c('Type','CG'),drop=FALSE],
         show_rownames=FALSE,show_colnames=FALSE,
         clustering_method = 'ward.D',
         cluster_cols = FALSE,
         annotation_colors=list(CellType=c('Jurkat'='#ff9289', 'K562'='#82b7ff', 'Mixed'='grey'),
                                Type=c('Aci'='grey','Hha Jurkat high'='#b3645eff','Hha K562 high'='#597daeff','significant'='purple')))
ggplot(to.plot,aes(x=CT,y=AMPL134269))+geom_violin()+theme_bw()
