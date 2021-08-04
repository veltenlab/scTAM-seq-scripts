library(caret)
#Train
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
train_control <- trainControl(method = "cv", number = 10)
 model <- train(CellType ~ .+CD3*CD33,
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
table(predict(model,train.data),train$CellType)

#Test
sample <- 'Sample3_default'
subsample <- FALSE
clust.file <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/cluster_assignment.tsv'))
dat.s3 <- read.table(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/',sample,'.barcode.cell.distribution.tsv'))
scrub <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores.csv'))
row.names(scrub) <- scrub$Barcode
ddoublet <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores_DoubletDetection.csv'))
row.names(ddoublet) <- row.names(dat.s3)
dcells <- read.csv(paste0('/users/mscherer/cluster/project/Methylome/analysis/missionbio/tapestri/',sample,'/tsv/doublet_scores_DoubletCells.csv'))
row.names(dcells) <- row.names(dat.s3)
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
                    #                   TotalProteinCounts=total.protein.counts,
                    CD3=prot.counts.cd3$raw,
                    CD33=prot.counts.cd33$raw)
res <- predict(model,test)
table(res,test$CellType)
