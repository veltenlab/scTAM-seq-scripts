library(ggplot2)
bm_cut <- 'GSM5935921_BM_HhaI'
bm_uncut <- 'GSM5935923_BM_undigested'
blood_cut <- 'GSM5935918_Blood_HhaI'
blood_uncut <- 'GSM5935920_Blood_undigested'
plot.path <- '~'
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_line(color='black', size=.25),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')
type_clear <- c("mut.only"="Mutation only",
                "NonHhaI"="No HhaI cutsite",
                "CpG.B.cell.diff"="Differential CpG Bcells",
                "CpG.MCL"="MCL-specific CpG",
                "CpG.MCL.multiple.cutsites"="MCL-specific CpG",
                "CpG.always.unmeth.B"="Always unmethylated CpG in Bcells",
                "CpG.B.cell.diff.and.MCL"="Differential CpG Bcells AND MCL-specific CpG",        
                "CpG.always.meth.B"="Always methylated CpG in Bcells",
                "CpG.imprinted.multiple.cutsites"="Imprinted CpG multiple",
                "CpG.imprinted"="Imprinted CpG")
colors_amplicons <- c("Mutation only"="#e5c494",
                      "Differential CpG Bcells"="#fc8d62",
                      "MCL-specific CpG"="#ffd92f",
                      "Differential CpG Bcells AND MCL-specific CpG"="#fc8d62",
                      "Always unmethylated CpG in Bcells"="#a6d854",
                      "Always methylated CpG in Bcells"="#66c2a5",
                      "Imprinted CpG"="#e78ac3",
                      "Imprinted CpG multiple"="#e78ac3",
                      "No HhaI cutsite"="#b3b3b3")
ampli_info <- read.table('../../misc/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
dat <- read.table(paste0("../../data/", bm_cut, "/tsv/", bm_cut, ".tsv.gz"), 
                      sep="\t", 
                      header=T)
rowinfo <- read.csv(paste0('../../misc/', bm_cut, '/tsv/rowinfo.csv'),
                    row.names = 1)
dropout_bm_cut <- apply(dat[row.names(rowinfo[which(rowinfo$DoubletDetectionLabel==0), ]), row.names(ampli_info)[grepl("NonHhaI", ampli_info$Type.of.amplicon)]], 2, function(x){
  sum(x==0)/length(x)
})
dat <- read.table(paste0("../../data/", bm_uncut, "/tsv/", bm_uncut, ".tsv.gz"), 
                  sep="\t", 
                  header=T)
rowinfo <- read.csv(paste0('../../misc/', bm_uncut, '/tsv/doublet_scores_DoubletDetection.csv'),
                    row.names = 2)
dropout_bm_uncut <- apply(dat[row.names(rowinfo[which(rowinfo$DoubletDetectionLabel==0), ]), row.names(ampli_info)[grepl("NonHhaI", ampli_info$Type.of.amplicon)]], 2, function(x){
  sum(x==0)/length(x)
})
dat <- read.table(paste0("../../data/", blood_cut, "/tsv/", blood_cut, ".tsv.gz"), 
                  sep="\t", 
                  header=T)
rowinfo <- read.csv(paste0('../../misc/', blood_cut, '/tsv/rowinfo.csv'),
                    row.names = 1)
dropout_blood_cut <- apply(dat[row.names(rowinfo[which(rowinfo$DoubletDetectionLabel==0), ]), row.names(ampli_info)[grepl("NonHhaI", ampli_info$Type.of.amplicon)]], 2, function(x){
  sum(x==0)/length(x)
})
dat <- read.table(paste0("../../data/", blood_uncut, "/tsv/", blood_uncut, ".tsv.gz"), 
                  sep="\t", 
                  header=T)
rowinfo <- read.csv(paste0('../../misc/', blood_uncut, '/tsv/rowinfo.csv'),
                    row.names = 1)
dropout_blood_uncut <- apply(dat[row.names(rowinfo[which(rowinfo$DoubletDetectionLabel==0), ]), row.names(ampli_info)[grepl("NonHhaI", ampli_info$Type.of.amplicon)]], 2, function(x){
  sum(x==0)/length(x)
})
to_plot <- data.frame(Blood_Undigested=dropout_blood_uncut,
                      Blood_Digested=dropout_blood_cut,
                      BM_Undigested=dropout_bm_uncut,
                      BM_Digested=dropout_bm_cut)
to_plot <- reshape2::melt(to_plot)
plot <- ggplot(to_plot, aes(x=variable, y=value))+geom_boxplot(color='black', fill='gray80', size=.25, outlier.size=.25)+plot_theme+xlab('Sample')+ylab('FNR')
ggsave(file.path(plot.path, 'boxplot_non_Hha_all.pdf'),
       width=75,
       height=37,
       unit='mm')
