############## plot_bottleneck_layer.R ##############
#' Using this script, one can visualize the bottleneck layer of the methylation autoencoder in a two dimensional UMAP or PCA
#' plot. Additional inputs to the script are the assignment of cells to cell types and the input read count file.

library(ggfortify)
sample <- 'BCells_Sample7_70_percent_good_performance'
bottleneck.layer <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/methylation_autoencoder/bottleneck.csv'),
                             header=TRUE)
bottleneck.layer <- bottleneck.layer[,-1]
input <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'),
                    sep='\t',
                    header=TRUE)
meth.data <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/methylation_autoencoder/mixture_prob_dca.csv'),
                      header = TRUE,
                      row.names = 1)
mis.amplis <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/methylation_autoencoder/missing_amplicons.csv'),
                       header = TRUE,
                       row.names = 1)
input <- input[,!(colnames(input)%in%as.character(mis.amplis[,1]))]
ampli.info <- read.table('/users/mscherer/cluster/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.tsv')
colnames(meth.data) <- colnames(input)
row.names(bottleneck.layer) <- row.names(meth.data) <- row.names(input)
pca.obj <- prcomp(bottleneck.layer)
autoplot(pca.obj)

library(M3C)
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
res <- tsne(t(bottleneck.layer),dotsize=1)
to.plot <- data.frame(res$data,AMPL130645=meth.data[,'AMPL130645'])
plot <- ggplot(to.plot,aes(x=X1,y=X2,color=AMPL130645))+geom_point()+
  ylab('tSNE2')+xlab('tSNE1')+theme+scale_color_gradientn(colors=rev(viridis(15)),
                                                              limits=c(0,1))

ggsave(paste0("/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/",sample,"/plots/bottleneck_layer_tSNE.pdf"),plot)

