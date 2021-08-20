library(caret)
#Train
subsample <- FALSE
sample <- 'Sample4_70_percent'
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
dat.s3 <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'))
scrub <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores.csv'))
row.names(scrub) <- scrub$Barcode
ddoublet <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores_DoubletDetection.csv'))
row.names(ddoublet) <- row.names(dat.s3)
dcells <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores_DoubletCells.csv'))
row.names(dcells) <- row.names(dat.s3)
dsolo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores_solo.csv'))
row.names(dsolo) <- row.names(dat.s3)
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
prot.counts <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample, '-protein-counts.tsv'),sep='\t',header = TRUE)
shared.barcodes <- intersect(prot.counts$cell_barcode,row.names(dat.s3))
prot.counts.cd3 <- prot.counts[prot.counts$ab_description=='CD3',]
prot.counts.cd33 <- prot.counts[prot.counts$ab_description=='CD33',]
row.names(prot.counts.cd33) <- prot.counts.cd33$cell_barcode
row.names(prot.counts.cd3) <- prot.counts.cd3$cell_barcode
prot.counts.cd33 <- prot.counts.cd33[shared.barcodes,]
prot.counts.cd3 <- prot.counts.cd3[shared.barcodes,]
prot.counts.cd3$raw[is.na(prot.counts.cd3$raw)] <- 0
prot.counts.cd33$raw[is.na(prot.counts.cd33$raw)] <- 0
clust.file <- clust.file[shared.barcodes,]
ddoublet <- ddoublet[shared.barcodes,]
ddoublet.scores <- ddoublet[,2]
ddoublet.scores[is.na(ddoublet.scores)] <- 0
dcells <- dcells[shared.barcodes,]
dcells.scores <- dcells[,2]
dsolo <- dsolo[shared.barcodes,]
dsolo.scores <- dsolo[,2]
dat.s3 <- dat.s3[shared.barcodes,]
scrub <- scrub[shared.barcodes,]
ampli.info <- read.csv('/users/mscherer/cluster/project/Methylome/infos/cell_lines/analysis_Agostina_Apr21/Summary_methylation_values_amplicons_K562_Jurkat_v2.csv')
dat.s3 <- dat.s3[as.character(ampli.info$AmpID_design_1898)]
#total.counts <- apply(dat.s3[,ampli.info$Type%in%"Aci"],1,sum)
total.counts <- apply(dat.s3,1,sum)
total.protein.counts <- prot.counts.cd3$raw+prot.counts.cd33$raw
train <- data.frame(CellType=as.factor(ifelse(clust.file$CellType=='Mixed',1,0)),
                   DoubletScore=scrub[row.names(clust.file),'DoubletScore'],
                   TotalCounts=total.counts[row.names(clust.file)],
                   DoubletDetectionScore=ddoublet.scores,
                   DoubletCellsScore=dcells.scores,
                   SoloScore=dsolo.scores,
#                   TotalProteinCounts=total.protein.counts,
                   CD3=prot.counts.cd3$raw,
                   CD33=prot.counts.cd33$raw)
#train <- na.omit(train)
if(subsample){
  sel.samples <- row.names(train[clust.file$CellType=='Mixed',])
  add.samples <- sample(row.names(train[clust.file$CellType!='Mixed',]),sum(clust.file$CellType=='Mixed'),replace=FALSE)
  train.data <- train[c(sel.samples,add.samples),]
}else{
  train$TotalCounts <- log10(train$TotalCounts+1)
  train$CD3 <- log10(train$CD3+1)
  train$CD33 <- log10(train$CD33+1)
#  train.data <- as.data.frame(apply(train[,-1],2,scale))
  train.data <- train
  train.data$CellType <- train$CellType
}
sample <- 'Sample3_80_percent'
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
dat.s3 <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'))
scrub <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores.csv'))
row.names(scrub) <- scrub$Barcode
ddoublet <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores_DoubletDetection.csv'))
row.names(ddoublet) <- row.names(dat.s3)
dcells <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores_DoubletCells.csv'))
row.names(dcells) <- row.names(dat.s3)
dsolo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores_solo.csv'))
row.names(dsolo) <- row.names(dat.s3)
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
prot.counts <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample, '-protein-counts.tsv'),sep='\t',header = TRUE)
shared.barcodes <- intersect(prot.counts$cell_barcode,row.names(dat.s3))
prot.counts.cd3 <- prot.counts[prot.counts$ab_description=='CD3',]
prot.counts.cd33 <- prot.counts[prot.counts$ab_description=='CD33',]
row.names(prot.counts.cd33) <- prot.counts.cd33$cell_barcode
row.names(prot.counts.cd3) <- prot.counts.cd3$cell_barcode
prot.counts.cd33 <- prot.counts.cd33[shared.barcodes,]
prot.counts.cd3 <- prot.counts.cd3[shared.barcodes,]
prot.counts.cd3$raw[is.na(prot.counts.cd3$raw)] <- 0
prot.counts.cd33$raw[is.na(prot.counts.cd33$raw)] <- 0
clust.file <- clust.file[shared.barcodes,]
ddoublet <- ddoublet[shared.barcodes,]
ddoublet.scores <- ddoublet[,2]
ddoublet.scores[is.na(ddoublet.scores)] <- 0
dcells <- dcells[shared.barcodes,]
dcells.scores <- dcells[,2]
dsolo <- dsolo[shared.barcodes,]
dsolo.scores <- dsolo[,2]
dat.s3 <- dat.s3[shared.barcodes,]
scrub <- scrub[shared.barcodes,]
ampli.info <- read.csv('/users/mscherer/cluster/project/Methylome/infos/cell_lines/analysis_Agostina_Apr21/Summary_methylation_values_amplicons_K562_Jurkat_v2.csv')
dat.s3 <- dat.s3[as.character(ampli.info$AmpID_design_1898)]
#total.counts <- apply(dat.s3[,ampli.info$Type%in%"Aci"],1,sum)
total.counts <- apply(dat.s3,1,sum)
total.protein.counts <- prot.counts.cd3$raw+prot.counts.cd33$raw
train.second <- data.frame(CellType=as.factor(ifelse(clust.file$CellType=='Mixed',1,0)),
                    DoubletScore=scrub[row.names(clust.file),'DoubletScore'],
                    TotalCounts=total.counts[row.names(clust.file)],
                    DoubletDetectionScore=ddoublet.scores,
                    DoubletCellsScore=dcells.scores,
                    SoloScore=dsolo.scores,
                    #                   TotalProteinCounts=total.protein.counts,
                    CD3=prot.counts.cd3$raw,
                    CD33=prot.counts.cd33$raw)
#train <- na.omit(train)
if(subsample){
  sel.samples <- row.names(train.second[clust.file$CellType=='Mixed',])
  add.samples <- sample(row.names(train.second[clust.file$CellType!='Mixed',]),sum(clust.file$CellType=='Mixed'),replace=FALSE)
  train.data <- train.second[c(sel.samples,add.samples),]
}else{
  train.second$TotalCounts <- log10(train.second$TotalCounts+1)
  train.second$CD3 <- log10(train.second$CD3+1)
  train.second$CD33 <- log10(train.second$CD33+1)
  #  train.data <- as.data.frame(apply(train[,-1],2,scale))
  train.data.second <- train.second
  train.data.second$CellType <- train.second$CellType
}
sample <- 'Sample2_70_percent'
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
dat.s3 <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'))
scrub <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores.csv'))
row.names(scrub) <- scrub$Barcode
ddoublet <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores_DoubletDetection.csv'))
row.names(ddoublet) <- row.names(dat.s3)
dcells <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores_DoubletCells.csv'))
row.names(dcells) <- row.names(dat.s3)
dsolo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores_solo.csv'))
row.names(dsolo) <- row.names(dat.s3)
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
prot.counts <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample, '-protein-counts.tsv'),sep='\t',header = TRUE)
shared.barcodes <- intersect(prot.counts$cell_barcode,row.names(dat.s3))
prot.counts.cd3 <- prot.counts[prot.counts$ab_description=='CD3',]
prot.counts.cd33 <- prot.counts[prot.counts$ab_description=='CD33',]
row.names(prot.counts.cd33) <- prot.counts.cd33$cell_barcode
row.names(prot.counts.cd3) <- prot.counts.cd3$cell_barcode
prot.counts.cd33 <- prot.counts.cd33[shared.barcodes,]
prot.counts.cd3 <- prot.counts.cd3[shared.barcodes,]
prot.counts.cd3$raw[is.na(prot.counts.cd3$raw)] <- 0
prot.counts.cd33$raw[is.na(prot.counts.cd33$raw)] <- 0
clust.file <- clust.file[shared.barcodes,]
ddoublet <- ddoublet[shared.barcodes,]
ddoublet.scores <- ddoublet[,2]
ddoublet.scores[is.na(ddoublet.scores)] <- 0
dcells <- dcells[shared.barcodes,]
dcells.scores <- dcells[,2]
dsolo <- dsolo[shared.barcodes,]
dsolo.scores <- dsolo[,2]
dat.s3 <- dat.s3[shared.barcodes,]
scrub <- scrub[shared.barcodes,]
ampli.info <- read.csv('/users/mscherer/cluster/project/Methylome/infos/cell_lines/analysis_Agostina_Apr21/Summary_methylation_values_amplicons_K562_Jurkat_v2.csv')
dat.s3 <- dat.s3[as.character(ampli.info$AmpID_design_1898)]
#total.counts <- apply(dat.s3[,ampli.info$Type%in%"Aci"],1,sum)
total.counts <- apply(dat.s3,1,sum)
total.protein.counts <- prot.counts.cd3$raw+prot.counts.cd33$raw
train.third <- data.frame(CellType=as.factor(ifelse(clust.file$CellType=='Mixed',1,0)),
                           DoubletScore=scrub[row.names(clust.file),'DoubletScore'],
                           TotalCounts=total.counts[row.names(clust.file)],
                           DoubletDetectionScore=ddoublet.scores,
                           DoubletCellsScore=dcells.scores,
                           SoloScore=dsolo.scores,
                           #                   TotalProteinCounts=total.protein.counts,
                           CD3=prot.counts.cd3$raw,
                           CD33=prot.counts.cd33$raw)
#train <- na.omit(train)
if(subsample){
  sel.samples <- row.names(train.third[clust.file$CellType=='Mixed',])
  add.samples <- sample(row.names(train.third[clust.file$CellType!='Mixed',]),sum(clust.file$CellType=='Mixed'),replace=FALSE)
  train.data <- train.third[c(sel.samples,add.samples),]
}else{
  train.third$TotalCounts <- log10(train.third$TotalCounts+1)
  train.third$CD3 <- log10(train.third$CD3+1)
  train.third$CD33 <- log10(train.third$CD33+1)
  #  train.data <- as.data.frame(apply(train[,-1],2,scale))
  train.data.third <- train.third
  train.data.third$CellType <- train.third$CellType
}
train_control <- trainControl(method = "cv", number = 10)
train.data <- as.data.frame(rbind(train.data,train.data.third))
train.data <- train.data[,c("CellType","SoloScore","DoubletDetectionScore","CD3","TotalCounts")]
 model <- train(CellType ~ .,
                data = train.data,
                trControl = train_control,
                method = "glm",
                family=binomial())
# model <- train(CellType ~ .,
#                data = train.data,
#                trControl = train_control,
#                method = "rf")
summary(model)
model$results
table(predict(model,train.data),train.data$CellType)

#Test
sample <- 'Sample5_80_percent'
subsample <- FALSE
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
dat.s3 <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'))
scrub <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores.csv'))
row.names(scrub) <- scrub$Barcode
ddoublet <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores_DoubletDetection.csv'))
row.names(ddoublet) <- row.names(dat.s3)
dcells <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores_DoubletCells.csv'))
row.names(dcells) <- row.names(dat.s3)
dsolo <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores_solo.csv'))
row.names(dsolo) <- row.names(dat.s3)
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
prot.counts <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample, '-protein-counts.tsv'),sep='\t',header = TRUE)
shared.barcodes <- intersect(prot.counts$cell_barcode,row.names(dat.s3))
prot.counts.cd3 <- prot.counts[prot.counts$ab_description=='CD3',]
prot.counts.cd33 <- prot.counts[prot.counts$ab_description=='CD33',]
row.names(prot.counts.cd33) <- prot.counts.cd33$cell_barcode
row.names(prot.counts.cd3) <- prot.counts.cd3$cell_barcode
prot.counts.cd33 <- prot.counts.cd33[shared.barcodes,]
prot.counts.cd3 <- prot.counts.cd3[shared.barcodes,]
prot.counts.cd3$raw[is.na(prot.counts.cd3$raw)] <- 0
prot.counts.cd33$raw[is.na(prot.counts.cd33$raw)] <- 0
clust.file <- clust.file[shared.barcodes,]
ddoublet <- ddoublet[shared.barcodes,]
ddoublet.scores <- ddoublet[,2]
ddoublet.scores[is.na(ddoublet.scores)] <- 0
dcells <- dcells[shared.barcodes,]
dcells.scores <- dcells[,2]
dsolo <- dsolo[shared.barcodes,]
dsolo.scores <- dsolo[,2]
dat.s3 <- dat.s3[shared.barcodes,]
scrub <- scrub[shared.barcodes,]
ampli.info <- read.csv('/users/mscherer/cluster/project/Methylome/infos/cell_lines/analysis_Agostina_Apr21/Summary_methylation_values_amplicons_K562_Jurkat_v2.csv')
dat.s3 <- dat.s3[as.character(ampli.info$AmpID_design_1898)]
#total.counts <- apply(dat.s3[,ampli.info$Type%in%"Aci"],1,sum)
total.counts <- apply(dat.s3,1,sum)
total.protein.counts <- prot.counts.cd3$raw+prot.counts.cd33$raw
test <- data.frame(CellType=as.factor(ifelse(clust.file$CellType=='Mixed',1,0)),
                    DoubletScore=scrub[row.names(clust.file),'DoubletScore'],
                    TotalCounts=total.counts[row.names(clust.file)],
                    DoubletDetectionScore=ddoublet.scores,
                    DoubletCellsScore=dcells.scores,
                    SoloScore=dsolo.scores,
                    #                   TotalProteinCounts=total.protein.counts,
                    CD3=prot.counts.cd3$raw,
                    CD33=prot.counts.cd33$raw)
test$TotalCounts <- log10(test$TotalCounts+1)
test$CD3 <- log10(test$CD3+1)
test$CD33 <- log10(test$CD33+1)
#  train.data <- as.data.frame(apply(train[,-1],2,scale))
res <- predict(model,test)
table(res,test$CellType)

# plotting
library(grid)
library(gridExtra)
library(ggplotify)
my_theme <- theme_bw()+theme(axis.text.x = element_text(size=1),
                             axis.text.y = element_text(size=1),
                             axis.title.x = element_text(size=3),
                             axis.title.y = element_text(size=3),
                             legend.position = 'none')
plots <- lapply(colnames(train)[-1], function(x){
  lapply(colnames(train)[-1],function(y){
    to.plot <- train[,c("CellType",x,y)]
    plot <- ggplot(to.plot,aes_string(x=x,y=y,color="CellType"))+geom_point(size=.1)+my_theme
    as.grob(plot)
  })})
pdf('/users/mscherer/cluster/project/Methylome/Rplot.pdf')
grid.arrange(grobs=unlist(plots,recursive = FALSE))
dev.off()

